use quick_xml::{Reader, Error, events::Event, name::QName};
use crate::metabolite::Metabolite;

use std::fs::File;
use std::io::{BufReader, Read};

fn parse_float_or_zero(s: &str) -> f64 {
    match s.parse::<f64>() {
        Ok(n) => n,
        Err(_) => 0.0,
    }
}

pub fn parse_xml(path: String) -> Result<Vec<Metabolite>, Error> {
    let file: File = File::open(path).unwrap();
    let mut buf_reader: BufReader<File> = BufReader::new(file);

    let mut buf: Vec<u8> = Vec::new();
    buf_reader.read_to_end(&mut buf).unwrap();
    let mut reader: Reader<&[u8]> = Reader::from_reader(&buf[..]);

    let mut metabolites: Vec<Metabolite> = Vec::new();

    while let Ok(event) = reader.read_event() {
        match event {
            Event::Start(ref e) if e.name() == QName(b"metabolite") => {
                let mut endo_exo: Vec<String> = Vec::new();

                let mut found_name: bool = false;
                let mut found_accession: bool = false;
                let mut found_smiles: bool = false;
                let mut found_mw: bool = false;
                let mut found_mf: bool = false;

                let mut accession: String = String::from("");
                let mut name: String = String::from("");
                let mut smiles: String = String::from("");
                let mut weight: String = String::from("");
                let mut formula: String = String::from("");

                loop {
                    match reader.read_event().unwrap() {
                        Event::Start(ref e) if e.name() == QName(b"accession") => {
                            if !found_accession {
                                accession = reader.read_text(QName(b"accession")).unwrap().to_string();
                                found_accession = true;
                            }
                        },
                        Event::Start(ref e) if e.name() == quick_xml::name::QName(b"name") => {
                            if !found_name {
                                name = reader.read_text(QName(b"name")).unwrap().to_string();
                                found_name = true;
                            }
                        },
                        Event::Start(ref e) if e.name() == quick_xml::name::QName(b"smiles") => {
                            if !found_smiles {
                                smiles = reader.read_text(QName(b"smiles")).unwrap().to_string();
                                found_smiles = true;
                            }
                        },
                        Event::Start(ref e) if e.name() == quick_xml::name::QName(b"monisotopic_molecular_weight") => {
                            if !found_mw {
                                weight = reader.read_text(QName(b"monisotopic_molecular_weight")).unwrap().to_string();
                                found_mw = true;
                            }
                        },
                        Event::Start(ref e) if e.name() == quick_xml::name::QName(b"chemical_formula") => {
                            if !found_mf {
                                formula = reader.read_text(QName(b"chemical_formula")).unwrap().to_string();
                                found_mf = true;
                            }
                        },
                        Event::Start(ref e) if e.name() == QName(b"ontology") => {
                            loop {
                                match reader.read_event().unwrap() {
                                    Event::Start(ref e) if e.name() == QName(b"descendants") => {
                                        loop {
                                            match reader.read_event().unwrap() {
                                                Event::Start(ref e) if e.name() == QName(b"descendant") => {
                                                    loop {
                                                        match reader.read_event().unwrap() {
                                                            Event::Start(ref e) if e.name() == QName(b"term") => {
                                                                let term: String = reader.read_text(QName(b"term")).unwrap().to_string();
                                                                if term == String::from("Endogenous") && !endo_exo.contains(&term) {
                                                                    endo_exo.push("Endogenous".to_string());
                                                                    
                                                                } else if term == String::from("Exogenous") && !endo_exo.contains(&term){
                                                                    endo_exo.push("Exogenous".to_string());
                                                                    
                                                                }    
                                                            },
                                                            _ => ()
                                                        }
                                                        break;
                                                    }
                                                },
                                                _ => ()
                                            }
                                            break;
                                        }
                                    },
                                    _ => ()
                                }
                                break;
                            }
                        }
                        Event::End(ref e) => {
                            if e.name() == QName(b"metabolite") {

                                let met: Metabolite = Metabolite {accession, 
                                                                  name, 
                                                                  formula,
                                                                  mz: parse_float_or_zero(&weight),
                                                                  smiles,
                                                                  endo_exo};
                                metabolites.push(met);
                                break;
                            }
                        }

                        Event::Eof => panic!("unexpected end of document"),
                        _ => (),
                    }
                }
            },
            Event::Eof => break, // reached the end of the document
            //Err(e) => panic!("error reading event: {:?}", e),
            _ => (), // ignore other events
        }
    }

    Ok(metabolites)
    
}

pub fn parse_only_accession(path: &String) -> Result<Vec<String>, Error>{
    let file: File = File::open(path).unwrap();
    let mut buf_reader: BufReader<File> = BufReader::new(file);
    //buf_reader.trim_text(true);

    let mut buf: Vec<u8> = Vec::new();
    buf_reader.read_to_end(&mut buf).unwrap();
    let mut reader: Reader<&[u8]> = Reader::from_reader(&buf[..]);

    let mut metabolites: Vec<String> = Vec::new();
        while let Ok(event) = reader.read_event() {
        match event {
            Event::Start(ref e) if e.name() == QName(b"metabolite") => {
                let mut found_accession: bool = false;
                let mut accession: String = String::from("");

                loop {
                    match reader.read_event().unwrap() {
                        Event::Start(ref e) if e.name() == QName(b"accession") => {
                            if !found_accession {
                                accession = reader.read_text(QName(b"accession")).unwrap().to_string();
                                found_accession = true;
                            }
                        },
                        Event::End(ref e) if e.name() == quick_xml::name::QName(b"metabolite") => {
                            // reached the end of the metabolite element

                            metabolites.push(accession);
                            break;
                        }
                        _ => (),
                    }
                }
            },
            Event::Eof => break, // reached the end of the document
            //Err(e) => panic!("error reading event: {:?}", e),
            _ => (), // ignore other events
        }
    }
    Ok(metabolites)
}

