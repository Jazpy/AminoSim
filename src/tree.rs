use crate::sequence::Sequence;
use crate::mutator::Mutator;

use std::collections::HashMap;

struct NNode {
    children: Vec<NNode>,
    id: Option<String>,
    branch_length: f64,
    sequence: Option<Sequence>
}

impl NNode {
    fn new_empty() -> NNode {
        NNode {
            children: Vec::<NNode>::new(),
            id: None,
            branch_length: 0.0,
            sequence: None
        }
    }

    fn set_id(&mut self, s: &str) {
        if s.len() > 0 {
            self.id = Some(String::from(s));
        } else {
            self.id = None;
        }
    }

    fn set_branch_length(&mut self, d: f64) {
        self.branch_length = d;
    }

    fn consume(&mut self, flag: u8, buf: &str) {
        match flag {
            1 => self.set_id(buf),
            2 => {
                let branch: f64 = match buf.parse() {
                    Ok(n)  => n,
                    Err(_) => panic!("Could not parse \"{}\" into branch!", buf)
                };

                self.set_branch_length(branch)
            }
            _ => assert!(false, "Invalid read flag = {}", flag)
        }
    }

    fn add_child(&mut self, c: NNode) {
        self.children.push(c);
    }

    #[allow(dead_code)]
    fn print(&self, indent_lvl: usize) {
        let indent: usize = indent_lvl * 1;

        println!("{:indent$}**********************", "", indent=indent);

        let mut stats = match &self.id {
            Some(s) => (*s).to_string(),
            None    => "NONE".to_string()
        };

        stats = format!("ID {} : Children {} : Branch Length {}",
            stats, self.children.len(), self.branch_length);

        println!("{:indent$}{}", "", stats, indent=indent);
        println!("{:indent$}children:", "", indent=indent);

        for child in &self.children {
            let child_id = match &child.id {
                Some(s) => (*s).to_string(),
                None    => "NONE".to_string()
            };

            println!("{:indent$}{}", "", child_id, indent=indent);
        }

        for child in &self.children {
            child.print(indent + 1);
        }
    }
}

pub struct NTree {
    root: Option<NNode>,
    size: usize,
    partition: usize,
    build_str: String
}

impl NTree {
    pub fn new(p: usize, s: String) -> NTree {
        NTree {
            root: None,
            size: 0,
            partition: p,
            build_str: s
        }
    }

    pub fn build_from_newick(&mut self) {
        assert!(self.root.is_none(), "Tree already built!");

        // Iterate over all chars, we'll use a stack to keep track of parent
        // nodes, and build the tree depth-first as we read new nodes.
        let mut iter  = self.build_str.chars();
        let mut stack = Vec::<NNode>::new();

        // String buffer to parse ids and branch lengths
        let mut buffer = "".to_string();
        let mut c_o    = iter.next();

        // Flag to keep track of what we're putting in the buffer
        //  1 - id, 2 - branch length
        let mut read_flag: u8 = 1;
        // Flag that indicates if we're done reading
        let mut break_bool    = false;

        // Node we're currently building, we'll start with the root
        let mut curr_node = NNode::new_empty();

        while let Some(c) = c_o {
            // If we find a new opening paren,
            // we're reading curr_node's grandchildren
            if c == '(' {
                stack.push(curr_node);
                curr_node = NNode::new_empty();
            // If we've read a comma or close paren, we've finished reading
            // a node, add it to its parent.
            } else if c == ',' || c == ')' {
                // Assert that there's a parent node to add to
                let stk_len = stack.len();
                assert!(stk_len > 0, "Empty tree building stack, does your
                    Newick tree have a single root node?");

                // Finish reading the current node
                curr_node.consume(read_flag, &mut buffer.trim());
                buffer.clear();
                read_flag = 1;

                // Finally, add the newly finished node to the its parent
                stack[stk_len - 1].add_child(curr_node);
                self.size += 1;

                // If we read a comma, keep reading at this level
                if c == ',' {
                    curr_node = NNode::new_empty();
                // Otherwise, pop up one level
                } else {
                    curr_node = stack.pop().unwrap();
                }
            // A colon delimits id and branch length
            } else if c == ':' {
                curr_node.consume(read_flag, &mut buffer.trim());
                buffer.clear();
                read_flag = 2;
            // Colon marks end of newick tree
            } else if c == ';'{
                curr_node.consume(read_flag, &mut buffer.trim());
                buffer.clear();
                break_bool = true;
            // Else, we're reading an id or branch length, put in buffer
            } else {
                buffer.push(c);
            }

            c_o = iter.next();

            // Check if we should break out of reading
            if break_bool {
                break;
            }
        }

        if c_o.is_some() {
            println!("Newick tree string included characters after
                ';' character. Ignoring...");
        }

        // Assert that the tree was paren balanced (no nodes left on stack)
        assert!(stack.len() == 0, "Unbalanced parens on Newick tree");
        self.root = Some(curr_node);
        self.size += 1;

        // Cleanup
        self.build_str = String::new();
    }

    pub fn dfs_evolve(&mut self, m: &dyn Mutator,
        h: &mut HashMap<String, Sequence>) {
        let mut curr_node = match &mut self.root {
            Some(root_node) => root_node,
            None            => panic!("Can't evolve an empty tree")
        };

        // Simple DFS, mutating ancestral as we advance through the tree
        let mut stack = Vec::<(&mut NNode, Option<&Sequence>)>::new();
        stack.push((curr_node, None));

        while !stack.is_empty() {
            let tuple = stack.pop().unwrap();
            curr_node = tuple.0;
            let parent_seq = tuple.1;

            // Build sequence for this node if it doesn't exist
            if parent_seq.is_some() {
                let mutated = m.mutate(parent_seq.unwrap(),
                    curr_node.branch_length);
                curr_node.sequence = Some(mutated);
            } else {
                assert!(curr_node.sequence.is_some(), "Can't evolve a tree
                    with no ancestral sequence");
            }

            // If no children, we reached a tip node and can add to result
            if curr_node.children.is_empty() {
                assert!(curr_node.id.is_some(), "Currently, only named tip
                    nodes are supported for evolution");
                h.insert((&curr_node.id.as_ref().unwrap()).to_string(),
                    curr_node.sequence.as_ref().unwrap().clone());
                continue
            }

            // Push all children with parent sequence (curr's sequence)
            for child in &mut curr_node.children {
                stack.push((child, curr_node.sequence.as_ref()));
            }
        }
    }

    pub fn create_ancestral(&mut self, m: &dyn Mutator) {
        let root = match &mut self.root {
            Some(r) => r,
            None    => panic!("Can't create ancestral for an empty tree")
        };

        root.sequence = Some(m.random(self.partition));
    }

    #[allow(dead_code)]
    pub fn print(&self) {
        match &self.root {
            Some(root_node) => root_node.print(0),
            None            => println!("Empty NTree")
        }
    }

    #[allow(dead_code)]
    pub fn get_size(&self) -> usize {
        self.size
    }

    #[allow(dead_code)]
    pub fn get_partition(&self) -> usize {
        self.partition
    }
}
