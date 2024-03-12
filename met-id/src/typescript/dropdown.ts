import { check_checkboxes, check_selected} from "./ms1/ms1_main";
import { invoke } from "@tauri-apps/api/tauri";


export function fill_dropdown(items: string[], elementID: string) {
    const dropdown = document.getElementById(elementID) as HTMLSelectElement;
    dropdown.innerHTML = "";
    for (let i = 0; i < items.length; i++) {
        const option = document.createElement("option");
        option.text = items[i];
        option.value = items[i];
        dropdown.add(option);
    }
}


const target_metabolome = {"HMDB (All)": ["Endogenous Metabolites", "Exogenous Metabolites", "Unspecified Metabolites"],
                           "HMDB (Brain)": ["Endogenous Metabolites", "Exogenous Metabolites", "Unspecified Metabolites"],
                           "HMDB (CSF)": ["Endogenous Metabolites", "Exogenous Metabolites", "Unspecified Metabolites"],
                           "Lipidmaps": []}

const metabolome_values = {"HMDB (All)": ["Endogenous", "Exogenous", "Unspecified"],
                           "HMDB (Brain)": ["Endogenous", "Exogenous", "Unspecified"],
                           "HMDB (CSF)": ["Endogenous", "Exogenous", "Unspecified"],
                           "Lipidmaps": []}


/*
const target_matrix = {"Positive mode": ["[M+H]+", "[M+K]+", "[M+Na]+"],
                       "Negative mode": ["[M-H]-", "[M+K-2H]-", "[M+Na-2H]-", "[M+Cl]-"], 
                       "FMP-10": ["Phenolic Hydroxyls", "Primary Amines"],
                       "AMPP": ["Aldehydes", "Carboxylic Acids"],
                       "Norharmane": ["[M+H]+", "[M+K]+", "[M+Na]+", "[M+Norharmane+H]+"],
                       }
*/


export async function new_tgt_matrix() {
    let new_trg_mtx: Record<string, string[]> = await invoke("matrix_dropdown_tauri", {});
    //console.log("new_trg_mtrx", new_trg_mtx)
    return new_trg_mtx
}




export const fill_under_dropdown = {
    
    template : 
            `<div class="checkbox-label-container">
            <p>{{value}}</p>
            <input type="checkbox" {{checked}} value={{id}} id={{id}}>
            </div>`,
    
    metabolites: function(dropdown_id: string, elementID: string) {
        const template = 
            `<div class="checkbox-label-container">
            <p>{{value}}</p>
            <input type="checkbox" {{checked}} value={{id}} id={{id}}>
            </div>`
    
        const myElement = document.getElementById(elementID);
        myElement!.innerHTML = "";
    
    
        const option_selected = (document.getElementById(dropdown_id)! as HTMLSelectElement).value;
    
    
        for (let i = 0; i < target_metabolome[option_selected].length; i++) {
            const renderedTemplate = template.replace("{{value}}", target_metabolome[option_selected][i]);
            const renderedTemplate3 = renderedTemplate.replaceAll("{{id}}", metabolome_values[option_selected][i]);
            const renderedTemplate2 = renderedTemplate3.replace("{{checked}}", 'checked="checked"');
            
            myElement!.innerHTML += renderedTemplate2;
        }
        for (let i = 0; i < target_metabolome[option_selected].length; i++) {
            document.getElementById(metabolome_values[option_selected][i])?.addEventListener("change", async () => {
                    
                let met_selected: string = check_selected("metabolome-dropdown");
                let matrix_selected: string = check_selected("matrix-dropdown");
    
                let met_type: string[] = check_checkboxes("ms1-metabolome-div1");
                let adducts: string[] = check_checkboxes("ms1-metabolome-div2");
                console.log(adducts);
                
                let count: number = await invoke("sql_counter_tauri", {met: met_selected, mat: matrix_selected, typ: met_type, adducts:adducts});
                let db_size = document.getElementById("db_size");
                
                if (db_size !== null) {
                    console.log("printing number")
                    db_size.textContent = count.toString();
                }
            });
        }
    },

    matrices: async function(dropdown_id: string, elementID: string) {
        let new_target_matrix = await new_tgt_matrix();
        console.log(new_target_matrix)
        const myElement = document.getElementById(elementID);
        myElement!.innerHTML = "";
        const option_selected = (document.getElementById(dropdown_id)! as HTMLSelectElement).value;
    
        for (let i = 0; i < new_target_matrix[option_selected].length; i++) {
            const renderedTemplate = this.template.replace("{{value}}", new_target_matrix[option_selected][i]);
            const renderedTemplate3 = renderedTemplate.replaceAll("{{id}}", `'${new_target_matrix[option_selected][i]}'`);
            let renderedTemplate2 = renderedTemplate3.replace("{{checked}}", 'checked="checked"');            
            myElement!.innerHTML += renderedTemplate2;
            
        }
        console.log(new_target_matrix[option_selected])

        let met_selected: string = check_selected("metabolome-dropdown");
        let matrix_selected: string = check_selected("matrix-dropdown");

        let met_type: string[] = check_checkboxes("ms1-metabolome-div1");
        let adducts: string[] = check_checkboxes("ms1-metabolome-div2");
        console.log(adducts);
        let count: number = await invoke("sql_counter_tauri", {met: met_selected, mat: matrix_selected, typ: met_type, adducts:adducts});
        let db_size = document.getElementById("db_size");
        
        if (db_size !== null) {
            console.log("printin number")
            db_size.textContent = count.toString();
        }
        for (let i = 0; i < new_target_matrix[option_selected].length; i++) {
            document.getElementById(new_target_matrix[option_selected][i])!.addEventListener("change", async () => {
                console.log("here")
                let met_selected: string = check_selected("metabolome-dropdown");
                let matrix_selected: string = check_selected("matrix-dropdown");
    
                let met_type: string[] = check_checkboxes("ms1-metabolome-div1");
                let adducts: string[] = check_checkboxes("ms1-metabolome-div2");
                console.log(adducts);
                let count: number = await invoke("sql_counter_tauri", {met: met_selected, mat: matrix_selected, typ: met_type, adducts:adducts});
                let db_size = document.getElementById("db_size");
                
                if (db_size !== null) {
                    console.log("printgin number")
                    db_size.textContent = count.toString();
                }
            });
        }
    }
}
