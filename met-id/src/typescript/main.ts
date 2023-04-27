import { invoke } from "@tauri-apps/api/tauri";
import { appWindow } from '@tauri-apps/api/window'
// Call the Tauri endpoint to open a file dialog and read the selected file
import { open } from '@tauri-apps/api/dialog';
// Open a selection dialog for image files
//import { draw } from './spectra';
import * as d3 from 'd3';
import { fill_dropdown, fill_options_under_dropdown } from "./dropdown";
import './ms1_table.ts';
import './ms1_sidebar.ts';



fill_dropdown(["HMDB (All)", "HMDB (Brain)", "HMDB (CSF)", "Lipidmaps"], "metabolome-dropdown");
fill_dropdown(["Positive mode", "Negative mode", "FMP-10", "AMPP", "Norharmane"], "matrix-dropdown");


fill_options_under_dropdown("metabolome", "metabolome-dropdown", "metabolome-checkbox-container")
fill_options_under_dropdown("matrix", "matrix-dropdown", "matrix-checkbox-container")


async function openFile() {
	const selected = await open({
	multiple: true,
	filters: [{
		name: 'MSMS_data',
		extensions: ['mzML']
	}]
	});
	if (Array.isArray(selected)) {
	// user selected multiple files
	} else if (selected === null) {
	// user cancelled the selection
	} else {
	// user selected a single file
	}
}


window.addEventListener("DOMContentLoaded", () => {
	document.getElementById("ms2-sidebar-open-file-button")?.addEventListener("click", () => openFile());

	document.getElementById('titlebar-minimize')?.addEventListener('click', () => appWindow.minimize());
	document.getElementById('titlebar-maximize')?.addEventListener('click', () => appWindow.toggleMaximize());
	document.getElementById('titlebar-close')?.addEventListener('click', () => appWindow.close());
});

function draw(div_index: number) {
	interface MyData {
		x: number;
		y: number;
		label: string;
	}
	
	const data: MyData[] = [
		{ x: 1, y: 2, label: 'Point 1' },
		{ x: 3, y: 2, label: 'Point 2' },
		{ x: 5, y: 6, label: 'Point 3' },
		// more data points here
	];




	const margin = { top: 20, right: 30, bottom: 30, left: 40 };
	const width = 900 - margin.left - margin.right;
	const height = 300 - margin.top - margin.bottom;

	var svg = d3.select("#results-card-spectrogram-div" + div_index as unknown as string)
	.append("svg")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
	.append("g")
		.attr("transform",
			"translate(" + margin.left + "," + margin.top + ")");

	const xScale = d3.scaleLinear()
			.domain([d3.min(data, d => d.x)!, d3.max(data, d => d.x)!])
			.range([0, width]);
/*
	let xAxis = svg.append("g")
		.attr("stroke", "white")
		.attr("stroke-width", 1)
		.attr("transform", "translate(0," + height + ")")
		.call(d3.axisBottom(xScale));
*/
	const yScale = d3.scaleLinear()
			.domain([d3.max(data, d => d.y)!, 0])
			.range([0, height]);
/*
	let yAxis = svg.append("g")
	.attr("stroke", "white")
	.attr("stroke-width", 1)
	.call(d3.axisLeft(yScale));
*/

	const line = d3.line<MyData>()
	.x(d => xScale(d.x))
	.y(d => yScale(d.y));

	svg.append('path')
	.datum(data)
	.attr('class', 'line')
	.attr('fill', 'none')
	.attr('stroke', 'white')
	.attr('stroke-width', 2.5)
	.attr('d', line);
}


const res_div = document.querySelector("#results-body") as HTMLDivElement;
const res_div_class = "results-card-spectrogram-div";

interface All_spectra_database {
	name: string[],
	//matrix: string[],
	smiles: string[],
	cid: string[],
	fragment: string[]
}
const all_spectra: All_spectra_database = await invoke("load_msms");

//const length_ = all_spectra[0].length;

function renderItems() {
	const batchSize = 3;
	let i = 0;

	function appendBatch() {
		const template = document.createElement("template");
		for (let j = 0; j < batchSize && i < 4; j++) {
		template.innerHTML = `
		<div class="results-card">
			<div class="results-card-textbox">
				<div class="results-card-name">
				<label class="results-card-name">${all_spectra[0][i]}</label>
				</div>
				<div class="results-card-dash">
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>CID</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>${all_spectra[2][i]}</p>
					</div>
				</div>
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>Matrix</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>FMP-10</p>
					</div>
				</div>
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>Adduct</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>${all_spectra[3][i]}</p>
					</div>
				</div>
				<div class="ms2-dash-box">
					<div class="ms2-dashbox-header">
					<p>Parent ion m/z</p>
					</div>
					<div class="ms2-dashbox-text">
					<p>400.000000</p>
					</div>
				</div>
				</div>
				<div class="molecule" id="molecule-canvas">
				<img data-smiles=${all_spectra[1][i]}
						data-smiles-options=" {'width': 300, 'height': 250 }"  
						data-smiles-theme='dark' />
				
				<script>
					SmiDrawer.apply({"themes": dark});
				</script>

				</div>
			</div>
			<div class="results-card-spectrogram-div" id="${res_div_class + i as unknown as string}"></div>

			</div>
			`

		const script = document.createElement('script');
		script.innerHTML = `
					SmiDrawer.apply();
					`;
		const templateContent = template.content;
		templateContent.appendChild(script);

		
		res_div.appendChild(templateContent);
		draw(i);
		i++
		}
		
		if (i < 4){//length_) {
			setTimeout(appendBatch, 0);
		}
	}
	appendBatch();
}

renderItems();


var slider = document.getElementById("myRange") as HTMLInputElement;
var output = document.getElementById("demo") as HTMLElement;
output!.innerHTML = slider!.value + "ppm";

slider!.oninput = function() {
	output!.innerHTML = (this as HTMLInputElement).value + " ppm";
}


const slideInDiv = document.querySelector('.slide-in');
const close_matrix_div_button = document.getElementById("close-add-matrix-div");
const add_buttons: string[] = ["ms1-sidebar-add-button", "ms1-sidebar-add-button2", "ms1-sidebar-add-button3"]

add_buttons.forEach((str) => {
	const add_matrix_button = document.getElementById(str);

	add_matrix_button!.addEventListener('click', () => {
		if (slideInDiv!.classList.value === "slide-in") {
			slideInDiv!.classList.add('slide-in--active');
		} else if (slideInDiv!.classList.value === "slide-in slide-in--active"){
			slideInDiv!.classList.remove('slide-in--active');
		}
		
	  });
	  close_matrix_div_button!.addEventListener('click', () => {
		  slideInDiv!.classList.remove('slide-in--active');
		});
});










































const button = document.getElementById("tab-2");

// Simulate a click event on the button
button!.click();



