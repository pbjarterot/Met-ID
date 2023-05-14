

export function fill_dropdown(items: string[], elementID: string) {
    const dropdown = document.getElementById(elementID) as HTMLSelectElement;
    dropdown.innerHTML = "";
    for (let i = 0; i < items.length; i++) {
        const option = document.createElement("option");
        option.text = items[i];
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


const target_matrix = {"Positive mode": ["[M+H]+", "[M+K]+", "[M+Na]+"],
                       "Negative mode": ["[M-H]-", "[M+K-2H]-", "[M+Na-2H]-", "[M+Cl]-"], 
                       "FMP-10": ["Phenols", "Primary amines", "Triols", "Diols", "Hydroxyls"],
                       "AMPP": ["Aldehydes", "Carboxylic Acids"]
                       }



//console.log((document.getElementById("metabolome-dropdown")! as HTMLSelectElement).value)

//doesnt trigger
export function fill_options_under_dropdown(type: string, dropdown_id: string, elementID: string) {
    const template = 
        `<div class="checkbox-label-container">
        <p>{{value}}</p>
        <input type="checkbox" checked="checked" value={{id}}>
        </div>`

    const myElement = document.getElementById(elementID);
    myElement!.innerHTML = "";


    const option_selected = (document.getElementById(dropdown_id)! as HTMLSelectElement).value;
    if (type === "metabolome")  {
        for (let i = 0; i < target_metabolome[option_selected].length; i++) {
            const renderedTemplate = template.replace("{{value}}", target_metabolome[option_selected][i]);
            const renderedTemplate2 = renderedTemplate.replace("{{id}}", metabolome_values[option_selected][i]);
            myElement!.innerHTML += renderedTemplate2;
        }
    } else if (type === "matrix") {
        for (let i = 0; i < target_matrix[option_selected].length; i++) {
            const renderedTemplate = template.replace("{{value}}", target_matrix[option_selected][i]);
            const renderedTemplate2 = renderedTemplate.replace("{{id}}", target_matrix[option_selected][i]);
                myElement!.innerHTML += renderedTemplate2;
        }
    }
}

