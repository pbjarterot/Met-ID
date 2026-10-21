use std::collections::{HashMap, HashSet};

#[derive(Debug, PartialEq)]
pub struct EquationCoefficients {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

pub fn parse_equation(input: &str) -> Option<EquationCoefficients> {
    let values: Vec<&str> = input.split(',').map(|s| s.trim()).collect();

    if input.is_empty() {
        return Some(EquationCoefficients {
            a: 0.0,
            b: 0.0,
            c: 0.0,
        });
    }

    match values.len() {
        1 => {
            if let Ok(constant) = values[0].parse::<f64>() {
                Some(EquationCoefficients {
                    a: 0.0,
                    b: 0.0,
                    c: constant,
                })
            } else {
                None
            }
        }
        2 => {
            if let (Ok(a), Ok(b)) = (values[0].parse::<f64>(), values[1].parse::<f64>()) {
                Some(EquationCoefficients { a, b, c: 0.0 })
            } else {
                None
            }
        }
        3 => {
            if let (Ok(a), Ok(b), Ok(c)) = (
                values[0].parse::<f64>(),
                values[1].parse::<f64>(),
                values[2].parse::<f64>(),
            ) {
                Some(EquationCoefficients { a, b, c })
            } else {
                None
            }
        }
        _ => None,
    }
}

pub fn calculate_y(coefficients: &EquationCoefficients, x: &f64) -> Option<f64> {
    let EquationCoefficients { a, b, c } = coefficients;
    let ppm_error: f64 = a * x * x + b * x + c;
    Some(x - (ppm_error * (x / 1_000_000.0)))
}

fn msms_exist(msms_set: &HashSet<&str>, query: &str) -> &'static str {
    if msms_set.contains(query) {
        "available"
    } else {
        "not-available"
    }
}

fn update_hashmaps(data: &mut Vec<Vec<HashMap<&'static str, String>>>) {
    let mut occurrences: HashMap<(String, String), i32> = HashMap::new(); // Map of (foo, bar) to its count
    let mut foo_counts: HashMap<String, i32> = HashMap::new(); // Map of foo to its total count

    // Count occurrences of each (foo, bar) pair and total foo counts
    for vec in data.iter() {
        for hashmap in vec.iter() {
            if let (Some(foo), Some(bar)) = (hashmap.get("names"), hashmap.get("matrix")) {
                *occurrences.entry((foo.clone(), bar.clone())).or_insert(0) += 1;
                *foo_counts.entry(foo.clone()).or_insert(0) += 1;
            }
        }
    }

    // Update the "coverage" value for each hashmap
    for vec in data.iter_mut() {
        for hashmap in vec.iter_mut() {
            if let Some(foo) = hashmap.get("names").map(String::as_str) {
                if let Some(total_foo_count) = foo_counts.get(foo) {
                    let mut count = *total_foo_count;
                    if let Some(bar) = hashmap.get("matrix").map(String::as_str) {
                        if let Some(same_pair_count) =
                            occurrences.get(&(foo.to_string(), bar.to_string()))
                        {
                            count -= same_pair_count - 1;
                        }
                    }
                    hashmap.insert("coverage", count.to_string());
                }
            }
        }
    }
}

#[tauri::command]
pub fn calculate_adjusted_mass(
    masses: Vec<String>,
    mass_error: String,
) -> Vec<HashMap<String, String>> {
    let mut input_masses: Vec<f64> = Vec::<f64>::new();
    let mut latent_variables: Vec<String> = Vec::<String>::new();
    for i in masses {
        match i.parse::<f64>() {
            Ok(num) => input_masses.push(num),
            Err(_e) => {
                println!("Treating {:?} as a latent variable", i);
                latent_variables.push(i);
            }
        }
    }

    let mut res: Vec<HashMap<String, String>> = Vec::new();
    let coefficients: Option<EquationCoefficients> = parse_equation(&mass_error);
    for mass in input_masses {
        let adj_mass: f64;
        if let Some(ref coeff) = coefficients {
            if let Some(result) = calculate_y(&coeff, &mass) {
                adj_mass = result;
            } else {
                println!("invalid_input");
                adj_mass = mass;
            }
        } else {
            println!("Invalid input");
            adj_mass = mass;
        }

        let mut res_map: HashMap<String, String> = HashMap::new();
        res_map.insert("oMass".to_string(), mass.to_string());
        res_map.insert("aMass".to_string(), adj_mass.to_string());
        res.push(res_map);
    }
    res
}

fn empty_match(mass: f64, adj_mass: f64, ccs: f64) -> HashMap<&'static str, String> {
    let mut res_map = HashMap::new();
    for key in [
        "tMass", "dMass", "dPPM", "names", "matrix", "accession",
        "smiles", "msms", "formula", "possible_derivs", 
    ] {
        res_map.insert(key, "".to_string());
    }
    res_map.insert("oMass", mass.to_string());
    res_map.insert("aMass", adj_mass.to_string());
    res_map.insert("ccs", ccs.to_string());
    res_map
}

// ---------- helpers: q-values from PEPs ----------
fn compute_qvals_exist(peps: &[(usize, f64)]) -> HashMap<usize, f64> {
    let mut items = peps.to_vec();
    items.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut qvals_raw: Vec<(usize, f64)> = Vec::with_capacity(items.len());
    let mut cum = 0.0;
    for (k, (idx, pep)) in items.iter().enumerate() {
        cum += pep;
        let q = cum / ((k + 1) as f64);
        qvals_raw.push((*idx, q));
    }
    // monotone from the bottom
    let mut qvals = qvals_raw.clone();
    let mut min_from_bottom = f64::INFINITY;
    for i in (0..qvals.len()).rev() {
        if qvals[i].1 < min_from_bottom { min_from_bottom = qvals[i].1; }
        qvals[i].1 = min_from_bottom;
    }
    qvals.into_iter().collect()
}

fn compute_qvals_top_hits(peps: &[(usize, usize, f64)])
    -> HashMap<(usize, usize), f64>
{
    let mut items = peps.to_vec();
    items.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

    let mut qvals_raw: Vec<((usize, usize), f64)> = Vec::with_capacity(items.len());
    let mut cum = 0.0;
    for (k, (pi, ri, pep)) in items.iter().enumerate() {
        cum += pep;
        let q = cum / ((k + 1) as f64);
        qvals_raw.push(((*pi, *ri), q));
    }
    let mut qvals = qvals_raw.clone();
    let mut min_from_bottom = f64::INFINITY;
    for i in (0..qvals.len()).rev() {
        if qvals[i].1 < min_from_bottom { min_from_bottom = qvals[i].1; }
        qvals[i].1 = min_from_bottom;
    }
    qvals.into_iter().collect()
}

fn background_density_uniform(adj_mass: f64, ppm_window: f64, sigma_ccs_eff: f64, k_ccs: f64, _k_mz: f64) -> f64 {
    // Treat the background as approximately uniform over the search
    // rectangle: m/z in [adj_mass ± delta] and CCS in ±0.5*sigma_ccs.
    // This only needs to be *relative*; units cancel in the posterior.
    
    let delta_mz = adj_mass * (ppm_window / 2.0) / 1e6;
    let ccs_span = 2.0 * k_ccs * sigma_ccs_eff; // ±3σ range
    let area = (2.0 * delta_mz) * ccs_span;
    if area > 0.0 { 1.0 / area } else { 1.0 }
}

// ========== Bivariate normal likelihood for candidates ==========
#[derive(Clone, Copy)]
struct MVN2 {
    mu_x: f64,
    mu_y: f64,
    inv_cov: Mat2,
    norm_const: f64, // 1 / (2*pi*sqrt(|Σ|))
}

impl MVN2 {
    fn from_params(mu_x: f64, mu_y: f64, sigma_x: f64, sigma_y: f64, rho: f64) -> Self {
        // Σ = [[σx^2, ρσxσy], [ρσxσy, σy^2]]
        let sxx = sigma_x * sigma_x;
        let syy = sigma_y * sigma_y;
        let sxy = rho * sigma_x * sigma_y;
        let cov = Mat2 { a11: sxx, a12: sxy, a21: sxy, a22: syy };
        let det = cov.det().abs().max(1e-30);
        let inv = cov.inv().unwrap_or(Mat2 { a11: 1e12, a12: 0.0, a21: 0.0, a22: 1e12 });
        let norm_const = 1.0 / (2.0 * std::f64::consts::PI * det.sqrt());
        MVN2 { mu_x, mu_y, inv_cov: inv, norm_const }
    }

    #[inline]
    fn pdf(&self, x: f64, y: f64) -> f64 {
        let dx = x - self.mu_x;
        let dy = y - self.mu_y;
        let q = self.inv_cov.quad_form(dx, dy);
        self.norm_const * (-0.5 * q).exp()
    }
}

// ========== small 2x2 math helpers ==========
#[derive(Clone, Copy, Debug)]
struct Mat2 { a11: f64, a12: f64, a21: f64, a22: f64 }

impl Mat2 {
    #[inline] fn det(&self) -> f64 { self.a11 * self.a22 - self.a12 * self.a21 }
    #[inline] fn inv(&self) -> Option<Mat2> {
        let d = self.det();
        if d.abs() < 1e-30 { return None; }
        Some(Mat2 {
            a11:  self.a22 / d,
            a12: -self.a12 / d,
            a21: -self.a21 / d,
            a22:  self.a11 / d,
        })
    }
    // (dx,dy)^T * self * (dx,dy)
    #[inline] fn quad_form(&self, dx: f64, dy: f64) -> f64 {
        dx * (self.a11 * dx + self.a12 * dy) + dy * (self.a21 * dx + self.a22 * dy)
    }
}

/* 
#[inline]
fn mean_cov_2d(xs: &[(f64, f64)]) -> (f64, f64, Mat2) {
    let n = xs.len().max(1) as f64;
    let (mut mx, mut my) = (0.0, 0.0);
    for &(x, y) in xs { mx += x; my += y; }
    mx /= n; my /= n;
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for &(x, y) in xs {
        let dx = x - mx; let dy = y - my;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    // ML covariance (divide by n). For unbiased, divide by n-1 when n>1.
    let cov = Mat2 { a11: sxx / n, a12: sxy / n, a21: sxy / n, a22: syy / n };
    (mx, my, cov)
}
*/

// New helper: q-values for the NULL decision (accepting H0)
fn compute_qvals_null(peps_null: &[(usize, f64)]) -> HashMap<usize, f64> {
    // Input: (peak_idx, PEP_null = 1 - P(H0|d))
    let mut items = peps_null.to_vec();
    items.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut qvals_raw: Vec<(usize, f64)> = Vec::with_capacity(items.len());
    let mut cum = 0.0;
    for (k, (idx, pepn)) in items.iter().enumerate() {
        cum += pepn;
        let q = cum / ((k + 1) as f64);
        qvals_raw.push((*idx, q));
    }
    // monotone from the bottom
    let mut qvals = qvals_raw.clone();
    let mut min_from_bottom = f64::INFINITY;
    for i in (0..qvals.len()).rev() {
        if qvals[i].1 < min_from_bottom { min_from_bottom = qvals[i].1; }
        qvals[i].1 = min_from_bottom;
    }
    qvals.into_iter().collect()
}

pub fn mass_matcher(
    input_masses: Vec<f64>,
    input_ccs: Vec<f64>,
    db_masses: &[f64],
    db_ccs: &[f64],
    db_names: Vec<String>,
    db_mnames: Vec<String>,
    db_accessions: Vec<String>,
    db_smiles: Vec<String>,
    db_formulas: Vec<String>,
    db_pos_derivs: Vec<i32>,
    mass_error: String,
    mass_window: f64,
    msms_ids: &Vec<String>,
    contains_ccs: bool,
) -> Vec<Vec<HashMap<&'static str, String>>> {

    // ---- Likelihood Σ for candidates (tune/estimate from QC) ----
    let sigma_mz: f64  = 0.00025;
    let _sigma_ppm: f64 = mass_window;
    let sigma_ccs: f64 = 2.0;
    let k_mz: f64 = 3.0;
    let k_ccs: f64 = 3.0;
    let rho: f64 = 0.0; //set > 0 if mz/ccs residuals are correlated
    let prediction_rmse: f64 = 5.103;
    let sigma_ccs_eff: f64 = (sigma_ccs.powf(2.0) + prediction_rmse.powf(2.0)).sqrt();

    // Prior for "none of the above"
    let pi0: f64 = 0.5;

    let mut res: Vec<Vec<HashMap<&'static str, String>>> = Vec::new();
    //let window: f64 = mass_window.parse::<f64>().unwrap();
    let ppm_factor: f64 = mass_window / 2.0 / 1e6;
    let coefficients: Option<EquationCoefficients> = parse_equation(&mass_error);
    
    let msms_set: HashSet<&str> = msms_ids.iter().map(|s| s.as_str()).collect();

    //Collect PEPs for q-values
    let mut exist_peps: Vec<(usize, f64)> = Vec::new();
    let mut top_hit_peps: Vec<(usize, usize, f64)> = Vec::new();

    for (idx, mass) in input_masses.iter().enumerate() {
        let adj_mass = calculate_y(coefficients.as_ref().unwrap(), mass).unwrap_or(*mass);
        
        let delta = adj_mass * ppm_factor;
        let (low_mass, high_mass) = (adj_mass - delta, adj_mass + delta);
        let low_idx = db_masses.partition_point(|&m| m < low_mass);
        let high_idx = db_masses.partition_point(|&m| m < high_mass);

        let cand_indices: Vec<usize> = (low_idx..high_idx).collect();
        if cand_indices.is_empty() {
            exist_peps.push((idx, 1.0));

            let mut res_vec: Vec<HashMap<&str, String>> = Vec::new();
            let mut row: HashMap<&str, String> = empty_match(*mass, adj_mass, input_ccs[idx]);

            row.insert("exist_pep", "1.0".to_string());
            row.insert("null_pep", "0.0".to_string());

            res_vec.push(row);
            res.push(res_vec);
            continue;
        }

        let mut likes = Vec::with_capacity(cand_indices.len());

        for &j in &cand_indices {
            //let sigma_mz_j = sigma_ppm * db_masses[j] / 1e6;
            let mvn = MVN2::from_params(
                db_masses[j],
                db_ccs[j],
                sigma_mz,
                sigma_ccs,
                rho
            );
            if contains_ccs {
                likes.push(mvn.pdf(adj_mass, input_ccs[idx]));
            }
            
        }

        //Background density via KDE
        //let u0 = kde.density(adj_mass, input_ccs[idx]);
        let u0 = background_density_uniform(adj_mass, mass_window, sigma_ccs_eff, k_ccs, k_mz);

        // Priors among windowed candidates
        let n_win = cand_indices.len() as f64;
        let pi_cand = if n_win > 0.0 { (1.0 - pi0) / n_win } else { 0.0 };

        //Evidence, posteriors, existence PEP
        let denom: f64 = pi0 * u0 + pi_cand * likes.iter().sum::<f64>();
        let p0 = (pi0 * u0) / denom;
        exist_peps.push((idx, p0));

        // Build rows
        let mut res_vec: Vec<HashMap<&'static str, String>> = Vec::with_capacity(high_idx - low_idx);
        let mut best_p: f64 = -1.0;
        let mut best_row_idx: usize = 0;

        for (loc, &j) in cand_indices.iter().enumerate() {
            //let mut like: f64 = 0.0;
            let mut post: f64 = 0.0;
            let mut ef: f64 = 0.0;
            if contains_ccs {
                let like = likes[loc];
                post = (pi_cand * like) / denom;
                ef = (1.0 - post) / post;

                if post > best_p { best_p = post; best_row_idx = loc; }
            }
            let t_mass = db_masses[j];
            let d_mass = adj_mass - t_mass;
            let t_ccs = db_ccs[j];
            let d_ppm  = (d_mass / t_mass) * 1e6;
            let mut d_ccs  = 0.0;
            if contains_ccs {
                d_ccs = input_ccs[idx] - t_ccs;
            }

            let db_name = if db_names[j].is_empty() {
                "Unknown Name"
            } else {
                &db_names[j]
            };

            let mut res_map = HashMap::new();
            res_map.insert("oMass", mass.to_string());
            res_map.insert("aMass", adj_mass.to_string());
            if contains_ccs {
                res_map.insert("ccs", input_ccs[idx].to_string());
                res_map.insert("tCCS", t_ccs.to_string());
                res_map.insert("ccsPPM", d_ccs.to_string());
            }
            res_map.insert("tMass", t_mass.to_string());
            res_map.insert("dMass", d_mass.to_string());
            res_map.insert("dPPM", d_ppm.to_string());
            res_map.insert("names", db_name.to_string());
            res_map.insert("matrix", db_mnames[j].clone());
            res_map.insert("accession", db_accessions[j].clone());
            res_map.insert("smiles", db_smiles[j].clone());
            res_map.insert("msms", msms_exist(&msms_set, &db_accessions[j]).to_string());
            res_map.insert("formula", db_formulas[j].clone());
            res_map.insert("possible_derivs", db_pos_derivs[j].to_string());


            //Bayes
            if contains_ccs {
                res_map.insert("posterior", post.to_string());
                res_map.insert("error_factor", ef.to_string());

                if p0 > best_p {
                    res_map.insert("exist_pep", p0.to_string());
                    res_map.insert("exist_q", "".to_string());
                } else {
                    let pep_top = 1.0 - best_p;
                    top_hit_peps.push((idx, best_row_idx, pep_top));
                }
            }
            

            res_vec.push(res_map);    
        }
        
        if res_vec.is_empty() {
            res_vec.push(empty_match(*mass, adj_mass, input_ccs[idx]));
        }
        res.push(res_vec);
    }
    // after building exist_peps = Vec<(peak_idx, p0)>
    let null_peps: Vec<(usize, f64)> = exist_peps.iter()
        .map(|(i, p0)| (*i, 1.0 - *p0))
        .collect();
    //compute overall q-values
    let exist_qvals: HashMap<usize, f64> = compute_qvals_exist(&exist_peps);
    let null_qvals: HashMap<usize, f64> = compute_qvals_null(&null_peps);
    let top_qvals: HashMap<(usize, usize), f64> = compute_qvals_top_hits(&top_hit_peps);

    // ---- Annotate rows with PEP/q-values ----
    for (peak_idx, rows) in res.iter_mut().enumerate() {
        if let Some(&q_exist) = exist_qvals.get(&peak_idx) {
            if let Some(&(_, p_exist)) = exist_peps.iter().find(|(i, _)| *i == peak_idx) {
                let q_null = *null_qvals.get(&peak_idx).unwrap_or(&1.0);
                let pep_null = 1.0 - p_exist;
                for row in rows.iter_mut() {
                    row.insert("exist_pep", p_exist.to_string()); // P(H0|d)
                    row.insert("exist_q", q_exist.to_string());

                    row.insert("null_pep", pep_null.to_string()); // 1 - P(H0|d)
                    row.insert("h0_q", q_null.to_string());
                }
            }
        }

        if let Some((&(pi, ri), &q_top)) = top_qvals.iter().find(|((pi, _ri), _)| *pi == peak_idx) {
            if pi == peak_idx && ri < rows.len() {
                if let Some(post_s) = rows[ri].get("posterior") {
                    if let Ok(post) = post_s.parse::<f64>() {
                        let pep = 1.0 - post;
                        rows[ri].insert("pep".into(), pep.to_string());
                        rows[ri].insert("q_value".into(), q_top.to_string());
                        rows[ri].insert("is_top".into(), "true".into());
                    }
                }
            }
        }
    }
    
    update_hashmaps(&mut res);
    /* 
    for i in &res {
        for j in i {
            println!("name: {:?}, posterior: {:?}, error_factor: {:?}, exist_pep: {:?}, exist_q: {:?}, pep: {:?}, q_value: {:?}, null pep: {:?}, null_q: {:?}", 
            j.get("names").unwrap(), 
            j.get("posterior"),
            j.get("error_factor"),
            j.get("exist_pep"),
            j.get("exist_q"),
            j.get("pep"),
            j.get("q_value"),
            j.get("null_pep"),
            j.get("h0_q")
        );
        }
    }
    */
    
    res
}
