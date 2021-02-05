mod parsers;
mod tree;
mod sequence;
mod mutator;

use crate::sequence::Sequence;

use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use clap::{Arg, App};

use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::prelude::*;

fn main() {
    // Get app info
    let matches = App::new("AminoSim")
        .version("0.9.0")
        .author("Jazeps Medina-Tretmanis <jaz.medtre@gmail.com>")
        .about("Fast amino acid simulation from coalescent trees.")
        .arg(Arg::with_name("treefile")
                 .short("t")
                 .long("treefile")
                 .takes_value(true)
                 .required(true)
                 .help("File with input coalescent tree(s)"))
        .arg(Arg::with_name("outfile")
                 .short("o")
                 .long("outfile")
                 .takes_value(true)
                 .required(true)
                 .help("Output filename"))
        .arg(Arg::with_name("length")
                 .short("l")
                 .long("length")
                 .takes_value(true)
                 .help("Length of generated sequences"))
        .arg(Arg::with_name("partitions")
                 .short("p")
                 .long("partitions")
                 .takes_value(true)
                 .help("File with coalescent tree partitions"))
        .arg(Arg::with_name("scale")
                 .short("s")
                 .long("scale")
                 .takes_value(true)
                 .help("Branch scaling factor"))
        .arg(Arg::with_name("threads")
                 .long("threads")
                 .takes_value(true)
                 .help("Maximum number of threads to spawn"))
        .get_matches();

    // Get args
    let tree_file = matches.value_of("treefile").unwrap();
    let out_file  = matches.value_of("outfile").unwrap();

    let partition_fp: Option<&str> = matches.value_of("partitions");

    let mut threads: usize = 1;
    let threads_arg = matches.value_of("threads");
    if threads_arg.is_some() {
        threads = match threads_arg.unwrap().parse::<usize>() {
            Ok(t) => t,
            Err(_) => panic!("--threads argument is not a positive integer")
        }
    }

    let mut scale: f64 = 1.0;
    let scale_arg = matches.value_of("scale");
    if scale_arg.is_some() {
        scale = match scale_arg.unwrap().parse::<f64>() {
            Ok(s) => s,
            Err(_) => panic!("--scale argument is not a float")
        }
    }

    // Initialize multithreading env
    ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

    // Parse coalescent tree inputs
    let parse_res = match partition_fp {
        Some(p) => parsers::parse_newick_partitioned(tree_file, p),
        None    => panic!("--length arg not implemented yet! Try --partitions")
    };

    let mut tree_vec = match parse_res {
        Ok(t)  => t,
        Err(x) => panic!("Parse error: {}", x)
    };

    println!("Done parsing trees");

    // Create a mutator model
    let mut_model = mutator::HKY::new(0.25, 0.25, 0.25, 0.25,
        'A' as u8, 'G' as u8, 'C' as u8, 'T' as u8, 1.0, scale);

    // Create ancestral sequences
    println!("Building ancestrals...");
    tree_vec.par_iter_mut().for_each(|t| t.create_ancestral(&mut_model));

    // Evolve all trees
    println!("Mutating ancestrals...");
    let mut mutated_seqs =
        vec![HashMap::<String, Sequence>::new(); tree_vec.len()];
    tree_vec.par_iter_mut().zip(mutated_seqs.par_iter_mut()).for_each(
        |(t, h)| t.dfs_evolve(&mut_model, h));
    tree_vec.clear();

    // Assemble mutant partitions
    println!("Assembling mutants...");
    let mut assembled_seqs = HashMap::<String, String>::new();
    for h in mutated_seqs {
        for (k, v) in h {
            let k_o = assembled_seqs.get_mut(&k);
            // If id exists in assembled sequences, append it
            if k_o.is_some() {
                k_o.unwrap().push_str(v.to_string())
            // If we haven't touched this id, add a new pair
            } else {
                assembled_seqs.insert(k, String::from(v.to_string())); ()
            }
        }
    }

    // Print out our mutants
    println!("Writing sequences...");
    let mut out = OpenOptions::new()
        .write(true)
        .create(true)
        .open(out_file)
        .unwrap();

    for (k, v) in assembled_seqs {
        if let Err(e) = writeln!(out, "{} {}", k, v) {
            panic!("Couldn't write to file: {}", e);
        }
    }

    println!("All done!");
}
