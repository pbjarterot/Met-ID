use log::info;
use once_cell::sync::OnceCell;
use std::path::PathBuf;

pub static METABOLITE_BIN_PATH: OnceCell<PathBuf> = OnceCell::new();

pub fn get_metabolite_bin_path() -> Result<PathBuf, std::io::Error> {
  let path: &PathBuf = METABOLITE_BIN_PATH.get().expect("Metabolite bin path cannot be found");
  info!("{:?}", path);
  Ok(path.to_owned())
}