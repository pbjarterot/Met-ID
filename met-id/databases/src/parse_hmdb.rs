use quick_xml::events::Event;
use quick_xml::{Reader, Error};
use crate::metabolite::Metabolite;


fn parse_float_or_zero(s: &str) -> f64 {
    match s.parse::<f64>() {
        Ok(n) => n,
        Err(_) => 0.0,
    }
}


pub fn parse_xml(path: String) -> Result<Vec<Metabolite>, Error> {
    let mut reader = Reader::from_file(path)?;
    reader.trim_text(true);

    let mut buf = Vec::new();
    let mut metabolites = Vec::new();
    loop {
        match reader.read_event(&mut buf) {
            Ok(Event::Start(ref e)) if e.name() == b"metabolite" => {
                let mut found_name = false;
                let mut found_accession = false;
                let mut found_smiles = false;
                let mut found_mw = false;
                let mut found_mf = false;

                let mut accession = String::from("");
                let mut name = String::from("");
                let mut smiles = String::from("");
                let mut weight = String::from("");
                let mut formula = String::from("");

                loop {
                    match reader.read_event(&mut buf) {
                        Ok(Event::Start(ref e)) if e.name() == b"accession" => {
                            if !found_accession {
                                accession = reader.read_text(b"accession", &mut Vec::new()).unwrap();
                                found_accession = true;
                            }
                        },
                        Ok(Event::Start(ref e)) if e.name() == b"name" => {
                            if !found_name {
                                name = reader.read_text(b"name", &mut Vec::new()).unwrap();
                                found_name = true;
                            }
                        },
                        Ok(Event::Start(ref e)) if e.name() == b"smiles" => {
                            if !found_smiles {
                                smiles = reader.read_text(b"smiles", &mut Vec::new()).unwrap();
                                found_smiles = true;
                            }
                        },
                        Ok(Event::Start(ref e)) if e.name() == b"monisotopic_molecular_weight" => {
                            if !found_mw {
                                weight = reader.read_text(b"monisotopic_molecular_weight", &mut Vec::new()).unwrap();
                                found_mw = true;
                            }
                        },
                        Ok(Event::Start(ref e)) if e.name() == b"chemical_formula" => {
                            if !found_mf {
                                formula = reader.read_text(b"chemical_formula", &mut Vec::new()).unwrap();
                                found_mf = true;
                            }
                        },
                        Ok(Event::End(ref e)) if e.name() == b"metabolite" => {
                            // reached the end of the metabolite element

                            let met = Metabolite {accession, 
                                                              name, 
                                                              formula,
                                                              mz: parse_float_or_zero(&weight),
                                                              smiles};
                            metabolites.push(met);
                            break;
                        },
                        Ok(Event::Eof) => panic!("unexpected end of document"),
                        _ => (),
                    }
                    buf.clear();
                }
            },
            Ok(Event::Eof) => break, // reached the end of the document
            Err(e) => panic!("error reading event: {:?}", e),
            _ => (), // ignore other events
        }
        buf.clear();
    }

    Ok(metabolites)
    
}
