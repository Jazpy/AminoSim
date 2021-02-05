use crate::sequence::Sequence;

use ndarray::arr2;

use std::f64::consts::E;
use rand::distributions::{Uniform, Distribution};

pub trait Mutator {
    fn mutate(&self, s: &Sequence, v: f64) -> Sequence;
    fn random(&self, l: usize) -> Sequence;
}

pub struct HKY {
    nuc_frequencies: [f64; 4],
    bases: [u8; 4],
    kappa: f64,
    beta: f64,
    scale: f64
}

impl HKY {
    pub fn new(pa: f64, pg: f64, pc: f64, pt: f64,
        ba: u8, bg: u8, bc: u8, bt: u8, k: f64, s: f64) -> HKY {
        // Calculate beta
        let b: f64 = 1.0 /
                     (2.0 * (pa + pg) * (pc + pt) +
                      2.0 * k * ((pa * pg) + (pc * pt)));

        HKY {
            nuc_frequencies: [pa, pg, pc, pt],
            bases: [ba, bg, bc, bt],
            kappa: k,
            beta: b,
            scale: s
        }
    }
}

impl Mutator for HKY {
    fn mutate(&self, s: &Sequence, v: f64) -> Sequence {
        let pa = self.nuc_frequencies[0];
        let pg = self.nuc_frequencies[1];
        let pc = self.nuc_frequencies[2];
        let pt = self.nuc_frequencies[3];

        let b = self.beta;
        let k = self.kappa;
        let scaled_v = v * self.scale;

        // TODO Move as much as possible to constructor
        let ag_ts_c = pa + pg + (pc + pt) * E.powf(-b * scaled_v);
        let ag_ts_e = E.powf(-(1.0 + (pa + pg) * (k - 1.0)) * b * scaled_v);
        let ct_ts_c = pc + pt + (pa + pg) * E.powf(-b * scaled_v);
        let ct_ts_e = E.powf(-(1.0 + (pc + pt) * (k - 1.0)) * b * scaled_v);
        let tv_c    = 1.0 - E.powf(-b * scaled_v);

        // Calculate A mutations
        let paa: f64 = (pa * ag_ts_c + pg * ag_ts_e) / (pa + pg);
        let pag: f64 = (pg * ag_ts_c - pg * ag_ts_e) / (pa + pg);
        let pac: f64 =  pc * tv_c;
        let pat: f64 =  pt * tv_c;

        // Calculate C mutations
        let pcc: f64 = (pc * ct_ts_c + pt * ct_ts_e) / (pc + pt);
        let pct: f64 = (pt * ct_ts_c - pt * ct_ts_e) / (pc + pt);
        let pca: f64 =  pa * tv_c;
        let pcg: f64 =  pg * tv_c;

        // Calculate G mutations
        let pgg: f64 = (pg * ag_ts_c + pa * ag_ts_e) / (pa + pg);
        let pga: f64 = (pa * ag_ts_c - pa * ag_ts_e) / (pa + pg);
        let pgc: f64 =  pc * tv_c;
        let pgt: f64 =  pt * tv_c;

        // Calculate T mutations
        let ptt: f64 = (pt * ct_ts_c + pc * ct_ts_e) / (pc + pt);
        let ptc: f64 = (pc * ct_ts_c - pc * ct_ts_e) / (pc + pt);
        let pta: f64 =  pa * tv_c;
        let ptg: f64 =  pg * tv_c;

        // Build matrix
        let matrix = arr2(&[
            [paa, pag, pac, pat],
            [pga, pgg, pgc, pgt],
            [pca, pcg, pcc, pct],
            [pta, ptg, ptc, ptt]
        ]);

        // Start mutating
        let mut mutated = s.nucleotides.clone();
        let mut rng = rand::thread_rng();
        let generator = Uniform::from(0.0..1.0);

        for n in mutated.iter_mut() {
            let row = if *n == self.bases[0] { 0 }
                else if  *n == self.bases[1] { 1 }
                else if  *n == self.bases[2] { 2 }
                else if  *n == self.bases[3] { 3 }
                else { panic!("Unrecognized base {} in Sequence being
                    mutated", n) };

            // Weighted random choice from transition probabilities
            let mut r: f64 = generator.sample(&mut rng);
            let mut new_base: u8 = 0;
            for i in 0..4 {
                let f = matrix[[row, i]];

                if r < f {
                    new_base = self.bases[i];
                    break
                }

                r -= f;
            }

            // Assert there's a valid new base
            assert!(new_base != 0, "Something went terribly wrong in Mutator's
                transition choice");

            // Modify the sequence with the new base
            *n = new_base;
        }

        // Build a Sequence object from mutated vec and freqs
        let mut freq_table = Vec::<(u8, f64)>::new();
        freq_table.push((self.bases[0], self.nuc_frequencies[0]));
        freq_table.push((self.bases[1], self.nuc_frequencies[1]));
        freq_table.push((self.bases[2], self.nuc_frequencies[2]));
        freq_table.push((self.bases[3], self.nuc_frequencies[3]));

        Sequence::from_vec(mutated, &freq_table)
    }

    fn random(&self, l: usize) -> Sequence {
        let mut freq_table = Vec::<(u8, f64)>::new();
        freq_table.push((self.bases[0], self.nuc_frequencies[0]));
        freq_table.push((self.bases[1], self.nuc_frequencies[1]));
        freq_table.push((self.bases[2], self.nuc_frequencies[2]));
        freq_table.push((self.bases[3], self.nuc_frequencies[3]));

        Sequence::new(&freq_table, l)
    }
}
