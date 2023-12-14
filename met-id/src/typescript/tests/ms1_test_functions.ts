import { invoke } from "@tauri-apps/api";


export async function count_metabolites(met_selected: string, matrix_selected: string, met_type: string[], adducts: string[]) {
  let count: number = await invoke("sql_counter", {met: met_selected, mat: matrix_selected, typ: met_type, adducts:adducts});

  return count;
}