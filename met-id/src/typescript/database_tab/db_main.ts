import { invoke } from '@tauri-apps/api';
import { CancellationToken } from '../ms2/ms2_main';

let DbCancellationToken: CancellationToken = { isCancelled: false };

window.addEventListener("DOMContentLoaded", async () => {
  renderDBHTML();

  const searchbar_container = document.getElementById("db-search");
  searchbar_container!.addEventListener("input", () => {
    console.log("clicked the searchbar");

    searchbar_container!.classList.remove("reanimate", "animate");

    setTimeout(() => {
      if ((searchbar_container! as HTMLInputElement).value == "") {
        animateDown(searchbar_container!);
      } else if (!searchbar_container!.classList.contains("animate")){
        animateUp(searchbar_container!);
      }
    }, 10)
    generate_results_cards((searchbar_container! as HTMLInputElement).value);
  })
});


async function generate_results_cards(inputvalue: string) {
  console.log(inputvalue);

  // Access the container div
  let div = document.getElementById("db-search-results-body");
  div!.innerHTML = "";
  if (!div) return; // Early exit if the div is not found

  let namelist2: Record<string, number> = await invoke("db_ids_and_names_tauri", {inputvalue});
  console.log("Here is namelist2: ", namelist2)

  for (const key in namelist2) {
    if (DbCancellationToken.isCancelled) {
      console.log("ive been cancelled")
      return;
    }
    //console.log("metname:", key, "id:", namelist2[key]);
    //const fragment = document.createDocumentFragment();
    const tempDiv = document.createElement('div');
    const L = await generate_db_small_results_card2(key, Number(namelist2[key]), false);
    tempDiv.innerHTML = `${L}`;
    tempDiv.classList.add("db-small-results-card")
    tempDiv.id = "db-small-results-card-" + namelist2[key];
    div!.appendChild(tempDiv);
    tempDiv!.addEventListener("click", async () => {
      //console.log("Ive been clicked");
      if (tempDiv.classList.contains("db-open")) {
        const L = await generate_db_small_results_card2(key, Number(namelist2[key]), false);
        tempDiv.innerHTML = `${L}`;
        tempDiv.classList.remove("db-open");
      } else {
        const L = await generate_db_small_results_card2(key, Number(namelist2[key]), true);
        tempDiv.innerHTML = `${L}`;
        tempDiv.classList.add("db-open");
      }
      let imdiv = document.getElementById("mol-im-container");

      // Generate and append the image asynchronously
      if (imdiv) {
        let image_area = await renderImageAsync();
        //console.log(image_area, imdiv);
        imdiv.appendChild(image_area);
      }

    })
  };

  // Locate the image container in the newly updated innerHTML
  let imdiv = document.getElementById("mol-im-container");
  //console.log(imdiv);

  // Generate and append the image asynchronously
  if (imdiv) {
    console.log("trying to make image")
    let image_area = await renderImageAsync();
    console.log(image_area, imdiv);
    imdiv.appendChild(image_area);
  }
}

export function renderDBHTML() {
  console.log("Hello");
  let db_search_div = document.getElementById("db-searchbar-div");

  db_search_div!.innerHTML = `
                <div class="db-content" id="db-content">
                  <div class='db-form-inline'>
                    <div class="db-input-group">
                      <input type='text' id='db-search' class="db-form-control db-search-form" placeholder="Search for a metabolite here">
                      <div class="db-search-results-body" id="db-search-results-body">
                      </div>
                    </div>
                  </div>
                </div>
                `
  
}

export async function generate_db_small_results_card2(name: string, index: number, open_: boolean) {
  //console.log("small_res2")

  let template = `<div class="db-small-results-card-top">
                    	<div class="db-small-results-card-textbox">
                      	<div class="db-small-results-card-name">
                        	<label class="db-small-results-card-name">${name}</label>
                    	  </div>
                  	  </div>
                	  </div>
                    <div class="db-molecule-information" id="db-molecule-information-${index}"></div>`
    if (open_) {
      console.log("index:", index);
      let dataHtml = await add_data_to_small_results_card(index);
      template = template.replace(`<div class="db-molecule-information" id="db-molecule-information-${index}"></div>`, 
                                  `<div class="db-molecule-information" id="db-molecule-information-${index}">${dataHtml}</div>`);
    }
    return template;
}

async function add_data_to_small_results_card(index: number){
  //console.log(index)
  let [name, formula, smiles, identifier, adductmap]: [string, string, string, string, Record<string, Record<string, number>>] = await invoke("db_data_tauri", {index})
  //console.log(adductmap);
  let temp  = ` <div class="db-molecule-information-top">
                  <div class="db-molecule-information">
                    <div class="db-molecule-name">${name}</div>
                    <div class="db-molecule-info-smaller" id="db-molecule-formula">${formula}</div>
                    <div class="db-molecule-info-smaller" id="db-molecule-id">${identifier}</div>
                    <div class="db-molecule-matrices">`
  
  for (const key in adductmap) {
    if (adductmap.hasOwnProperty(key)) {
      temp += adduct_div(key, adductmap[key])
    }
  }
  //adduct_div("AMPP", adductmap['AMPP'])}
  console.log(smiles)
  temp +=  `</div>
            </div>
            <div class="mol-im-container" id="mol-im-container">
            <img data-smiles="${smiles}" data-smiles-options="{'width': 600, 'height': 500, 'padding': 0.0}" data-smiles-theme='dark' />
            </div>
          </div>`            
  return temp
}

function adduct_div(name: string, a: Record<string, number>) {
  console.log(a)
  a = sortRecordByValues(a);
  let template = `<div class="db-adduct-div">	
                        <div class="db-adduct-div-name">
                          ${name}
                        </div>
                        <div class="db-adduct-div-matrix-div">`
  for (const key in a) {
    if (a.hasOwnProperty(key)) {
      template += adduct_card(key, a[key])
    }
  }
  template +=  `</div>
                      </div>`
  return template
}

function adduct_card(adduct: string, mass: number) {
	let template = `
	<div class="adduct-card">
		<div class="adduct-card-container">
			<h4><b>${adduct}</b></h4>
			<p>${roundToSignificantFigures(mass, 10)}</p>
		</div>
	</div>
	`
	return template
}
function roundToSignificantFigures(num: number, sigFigs: number): number {
  if (num === 0) return 0; // Early exit if the number is zero

  // Calculate the scale factor to move decimal point after sigFigs digits
  const scale = Math.pow(10, sigFigs - Math.ceil(Math.log10(Math.abs(num))));

  // Scale the number, round it, and then scale back
  return Math.round(num * scale) / scale;
}

function sortRecordByValues(record: Record<string, number>): Record<string, number> {
  // Convert record to an array of [key, value] tuples
  const entries = Object.entries(record);

  // Sort entries based on the value
  entries.sort((a, b) => a[1] - b[1]);

  // Convert the sorted array back to a record
  return entries.reduce((sortedRecord, [key, value]) => {
      sortedRecord[key] = value;
      return sortedRecord;
  }, {} as Record<string, number>);
}

async function renderImageAsync(): Promise<HTMLElement> {

    const molecule = document.createElement("div");
    molecule.className = "ms1-table-molecule";

    const imcontainer = document.createElement("div");
    imcontainer.className = "db-im-container";

    //imcontainer.innerHTML = `<img data-smiles=${smile_for_img} data-smiles-options="{'width': 600, 'height': 500, 'padding': 0.0}" data-smiles-theme='dark' />`;

    //const imcontainer2 = document.createElement("div");
    //imcontainer2.className = "db-im-container2";
    //imcontainer2.innerHTML = `<p class="ms1-molecule-name" onclick="window.open('https://hmdb.ca/metabolites/${hmdb_id}', '_blank')">${name_for_img}</p> `;

    const script = document.createElement("script");
    script.innerHTML = `SmiDrawer.apply();`;
    imcontainer.appendChild(script);
    molecule.appendChild(imcontainer);
    //molecule.appendChild(imcontainer2);
    return molecule;
}



function animateUp(element: HTMLElement): void {
  // Ensure any previous transition is cleared by removing the style property
  element.style.removeProperty('transform');
  // Trigger reflow to ensure the removal takes effect immediately
  void element.offsetWidth;

  // Apply transition and transform properties to move the element up
  element.style.transition = 'transform 1s';
  element.style.transform = 'translateY(-480px)';
  let a = document.getElementById("db-search-results-body") as HTMLDivElement;
  a.style.visibility = "visible";
  a.style.transition = "transform 1.3s";
  a.style.transform = "translateY(-480px)";
}




function animateDown(element: HTMLElement): void {
  // Ensure any previous transition is cleared by removing the style property
  element.style.removeProperty('transform');
  // Trigger reflow to ensure the removal takes effect immediately
  void element.offsetWidth;

  // Apply transition and transform properties to move the element back to its original position
  element.style.transition = 'transform 1s';
  element.style.transform = 'translateY(0)';

  let a = document.getElementById("db-search-results-body") as HTMLDivElement;
  a.style.visibility = "hidden";
  a.style.transition = 'transform 1s';
  a.style.transform = 'translateY(0)';
}
