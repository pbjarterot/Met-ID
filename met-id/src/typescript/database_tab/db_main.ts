import { invoke } from '@tauri-apps/api/core';


window.addEventListener("DOMContentLoaded", async () => {
  renderDBHTML();

  const searchbar_container = document.getElementById("db-search");
  searchbar_container!.addEventListener("input", () => {

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
  // Access the container div
  let div = document.getElementById("db-search-results-body");
  div!.innerHTML = "";
  if (!div) return; // Early exit if the div is not found

  let namelist2: [string, [string, number], number] = await invoke("db_ids_and_names_tauri", {inputvalue});

  for (let i = 0; i < namelist2.length; i++) {
    const tempDiv = document.createElement('div');
    const L = await generate_db_small_results_card2(namelist2[i][0], Number(namelist2[i][1][1]), namelist2[i][1][0], false);
    tempDiv.innerHTML = `${L}`;
    tempDiv.classList.add("db-small-results-card")
    tempDiv.id = "db-small-results-card-" + namelist2[i][1][1];
    div!.appendChild(tempDiv);

    function attachClickListener() {
        const childElement = document.getElementById('db-small-results-card-top-' + namelist2[i][1][1]);
        if (childElement) {
            childElement.addEventListener("click", async () => {
                if (tempDiv.classList.contains("db-open")) {
                    const L = await generate_db_small_results_card2(namelist2[i][0], Number(namelist2[i][1][1]), namelist2[i][1][0], false);
                    tempDiv.innerHTML = `${L}`;
                    tempDiv.classList.remove("db-open");
                } else {
                    const L = await generate_db_small_results_card2(namelist2[i][0], Number(namelist2[i][1][1]), namelist2[i][1][0], true);
                    tempDiv.innerHTML = `${L}`;
                    tempDiv.classList.add("db-open");
                }
                attachClickListener(); // Reattach the click listener to the new element

                let imdiv = document.getElementById("mol-im-container");
                // Generate and append the image asynchronously
                if (imdiv) {
                    let image_area = await renderImageAsync();
                    imdiv.appendChild(image_area);
                }
            });
        } else {
            console.error("Child element not found");
        }
    }
    attachClickListener();  // Attach the click listener for the first time
  };

  // Locate the image container in the newly updated innerHTML
  let imdiv = document.getElementById("mol-im-container");

  // Generate and append the image asynchronously
  if (imdiv) {
    let image_area = await renderImageAsync();
    imdiv.appendChild(image_area);
  }
}

export function renderDBHTML() {
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

export async function generate_db_small_results_card2(name: string, index: number, origin: string, open_: boolean) {

	let template = `<div class="db-small-results-card-top" id="db-small-results-card-top-${index}">
                    <div class="db-small-results-card-textbox">
                    	<div class="db-small-results-card-name">
                		<label class="db-small-results-card-name">${name}</label>
                    	</div>
					</div>
                	</div>
                    <div class="db-molecule-information" id="db-molecule-information-${index}"></div>`
    if (open_) {
    	let dataHtml = await add_data_to_small_results_card(index, origin);
    	template = template.replace(`<div class="db-molecule-information" id="db-molecule-information-${index}"></div>`, 
                                  `<div class="db-molecule-information" id="db-molecule-information-${index}">${dataHtml}</div>`);
    }
    return template;
}

export interface ParsedDBData {
  name: string;
  mz: string;
  db_accession: string;
  smiles: string;
  formula: string;
  map: Record<string, Record<string, number>>;
  functional_groups: Record<string, number>;
}

async function add_data_to_small_results_card(index: number, origin: string){
	let d: ParsedDBData = await invoke("db_data_tauri", {index, origin})
	console.log(d);
	if (origin == "lipids") {
		d.db_accession = "";
	}
	let temp  = `<div class="db-molecule-information-top">
				<div class="db-molecule-information">
				<div class="db-molecule-name">${d.name}</div>
				<div class="db-molecule-info-smaller" id="db-molecule-formula">${d.formula}</div>
				<div class="db-molecule-info-smaller" id="db-molecule-mass">${d.mz}</div>
				<div class="db-molecule-info-smaller" id="db-molecule-id">${d.db_accession}</div>
				<div class="db-molecule-matrices">`

	for (const key in d.map) {
		if (d.map.hasOwnProperty(key)) {
			temp += adduct_div(key, d.map[key])
		}
	}
	let functional_group_temp = ""

	for (const key in d.functional_groups) {
		if (d.functional_groups.hasOwnProperty(key) && key != "id" && d.functional_groups[key] != 0) {
			functional_group_temp += `<div>${key} = ${d.functional_groups[key]}</div>`
		}
	}
	console.log("funcional groups: ", functional_group_temp);


	temp +=  `</div>
			</div>
			<div class="mol-im-container" id="mol-im-container">
				<img data-smiles="${d.smiles}" data-smiles-options="{'width': 600, 'height': 500, 'padding': 0.0}" data-smiles-theme='dark' />
				${functional_group_temp}
			</div>
			

			</div>`            
	return temp
}

function adduct_div(name: string, a: Record<string, number>) {
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
