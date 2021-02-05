use crate::tree;

use rayon::prelude::*;

use std::fs::File;
use std::path::Path;
use std::io::{Result, Lines, BufReader, BufRead,
              stdout, Error, ErrorKind, Write};

fn read_lines<P>(filename: P) ->
    Result<Lines<BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

pub fn parse_newick_partitioned<P>(tree_fp: P, part_fp: P) ->
    Result<Vec::<tree::NTree>>
where P: AsRef<Path>, {
    // Iterators
    let mut tree_lines = read_lines(tree_fp)?;
    let mut part_lines = read_lines(part_fp)?;
    let iter = tree_lines.by_ref().zip(part_lines.by_ref());
    // Stats
    let mut line_counter: usize = 0;
    let mut part_counter: usize = 0;
    // Results
    let mut tree_vec = Vec::<tree::NTree>::new();

    for (tree_line_o, part_line_o) in iter {
        let tree_line = tree_line_o?;
        let part_line = part_line_o?;

        // First, try and parse the partition number
        let part: usize = match part_line.parse::<usize>() {
            Ok(n) => n,
            Err(_) => return Err(Error::new(ErrorKind::Other,
                format!("Could not parse partition '{}' into number",
                    part_line)))
        };

        part_counter += part;

        // Now that we have a partition length, create preliminary tree objs
        let tree_line = tree_line.trim();
        assert!(tree_line.ends_with(';'),
            "Incorrect Newick tree format, missing trailing ';'");

        let tree = tree::NTree::new(part, String::from(tree_line));
        tree_vec.push(tree);

        line_counter += 1;
        print!("\rDone reading {} trees and partitions", line_counter);
    }

    // Parse all trees in vector
    println!("\nParsing {} trees that cover {} bases...",
        line_counter, part_counter);
    stdout().flush()?;
    tree_vec.par_iter_mut().for_each(|t| t.build_from_newick());

    Ok(tree_vec)
}
