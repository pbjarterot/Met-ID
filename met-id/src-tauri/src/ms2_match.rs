

use rayon::prelude::*;

fn normalize_vector(input: &[i64]) -> Vec<f64> {
    let max_value = *input.iter().max().unwrap() as f64;
    input.iter().map(|&x| x as f64 / max_value).collect()
}

fn bin_spectrum(spectrum: &[(f64, f64)], bin_size: f64, threshold: f64) -> Vec<f64> {
    let mut binned_spectrum: Vec<f64> = Vec::new();

    if threshold != 0.0 {
        let max_mz: f64 = spectrum.iter().map(|(mz, _)| *mz).filter(|&mz| mz <= threshold).fold(f64::NAN, f64::max);
        let num_bins: usize = (max_mz / bin_size).ceil() as usize;
        binned_spectrum = vec![0.0; num_bins];

        for &(mz, intensity) in spectrum {
            if mz <= threshold {
                let index = (mz / bin_size) as usize;
                binned_spectrum[index] += intensity;
            }
        }
    } else {
        let max_mz: f64 = spectrum.iter().map(|(mz, _)| *mz).fold(f64::NAN, f64::max);
        let num_bins: usize = (max_mz / bin_size).ceil() as usize;
        binned_spectrum = vec![0.0; num_bins];
    
        for &(mz, intensity) in spectrum {
            let index = (mz / bin_size) as usize;
            binned_spectrum[index] += intensity;
        }
    }
    
    binned_spectrum
}

fn cosine_similarity(spectrum1: &[f64], spectrum2: &[f64]) -> f64 {
    let dot_product: f64 = spectrum1.iter().zip(spectrum2).map(|(a, b)| a * b).sum();
    let norm1: f64 = spectrum1.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
    let norm2: f64 = spectrum2.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
    dot_product / (norm1 * norm2)
}

pub fn ms2_matcher(query_spectrum: &(Vec<f64>, Vec<i64>), db_spectra: &[(Vec<f64>, Vec<i64>)], binsize: f64, threshold: f64) -> Vec<f64> {
    let normalized_query_intensity = normalize_vector(&query_spectrum.1);
    let spectrum1: Vec<(f64, f64)> = query_spectrum.0.iter().zip(normalized_query_intensity.iter()).map(|(&mz, &intensity)| (mz, intensity)).collect();
    let binned_spectrum1 = bin_spectrum(&spectrum1, binsize, threshold);

    db_spectra.par_iter().map(|spectrum| {
        let normalized_intensity = normalize_vector(&spectrum.1);
        let spectrum2: Vec<(f64, f64)> = spectrum.0.iter().zip(normalized_intensity.iter()).map(|(&mz, &intensity)| (mz, intensity)).collect();
        let binned_spectrum2 = bin_spectrum(&spectrum2, binsize, threshold);
        cosine_similarity(&binned_spectrum1, &binned_spectrum2)
    }).collect()
}