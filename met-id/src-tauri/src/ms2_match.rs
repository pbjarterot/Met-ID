

fn normalize_vector(input: &Vec<i64>) -> Vec<f64> {
    let &max_value = input.iter().max().unwrap();
    input.iter().map(|&x| x as f64 / max_value as f64).collect()
}

fn bin_spectrum(spectrum: &Vec<(f64, f64)>, bin_size: f64) -> Vec<f64> {
    let max_mz: f64 = spectrum.iter().map(|(mz, _)| *mz).fold(0.0/0.0, f64::max);
    let num_bins: usize = (max_mz / bin_size).ceil() as usize;
    let mut binned_spectrum: Vec<f64> = vec![0.0; num_bins];

    for (mz, intensity) in spectrum {
        let index: usize = (*mz / bin_size) as usize;
        binned_spectrum[index] += intensity;
    }
    binned_spectrum
}

fn cosine_similarity(spectrum1: &Vec<f64>, spectrum2: &Vec<f64>) -> f64 {
    let dot_product: f64 = spectrum1.iter().zip(spectrum2).map(|(a, b)| a * b).sum::<f64>();
    let norm1: f64 = spectrum1.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
    let norm2: f64 = spectrum2.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
    dot_product / (norm1 * norm2)
}

pub fn ms2_matcher(query_spectrum: (Vec<f64>, Vec<i64>), db_spectra: Vec<(Vec<f64>, Vec<i64>)>, binsize: f64) ->  Vec<f64> {
	let spectrum1: (Vec<f64>, Vec<f64>) = (query_spectrum.to_owned().0, normalize_vector(&query_spectrum.1));
	let spectrum1: Vec<(f64, f64)> = spectrum1.0.into_iter().zip(spectrum1.1.into_iter()).collect();

    //let bin_size: f64 = 0.01; //Daltons
    let binned_spectrum1: Vec<f64> = bin_spectrum(&spectrum1, binsize);

	let mut similarity_vec: Vec<f64> = Vec::new();

    for (index, spectrum) in db_spectra.iter().enumerate() {
		let spectrum2: Vec<(f64, f64)> = spectrum.to_owned().0.into_iter().zip(normalize_vector(&spectrum.1).to_owned().into_iter()).collect();
		
        let binned_spectrum: Vec<f64> = bin_spectrum(&spectrum2, binsize);
        let similarity: f64 = cosine_similarity(&binned_spectrum1, &binned_spectrum);
		similarity_vec.insert(index, similarity);
    }
	similarity_vec
}