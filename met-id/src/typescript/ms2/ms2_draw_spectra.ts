import * as d3 from 'd3';
import { invoke } from "@tauri-apps/api/tauri";
//import { zoom as d3Zoom, zoomIdentity as d3ZoomIdentity } from 'd3-zoom';

export interface SpectrumPoint {
    x: number;
    y: number;
}




export type SpectrumMap = {[key: string]: SpectrumPoint[]}


export async function draw(index: number, identifier: string, adduct: string, cid: string) {
    const fetchData = async (identifier: string, adduct: string, cid: string): Promise<{ keys: string[], values: SpectrumPoint[][] }> => {
        const result = await invoke<string>('get_msms_spectra_tauri', { identifier, adduct, cid });
        const resultMap = JSON.parse(result) as SpectrumMap;

        const keys: string[] = [];
        const values: SpectrumPoint[][] = [];

        for (const key in resultMap) {
            keys.push(key);
            values.push(resultMap[key]);
        }

        return { keys, values };
    };

	

    const colors = ["red", "green", "blue", "purple", "orange", "yellow"];
    const {keys: labels, values: dat} = await fetchData(identifier, adduct, cid);
		const shouldIncludeArray: boolean[] = new Array(labels.length).fill(true);

    const margin = { top: 40, right: 30, bottom: 30, left: 70 };
    const width = 1550 - margin.left - margin.right;
    const height = 420 - margin.top - margin.bottom;

    const svgContainer = d3.select("#ms2-small-results-card-plot-" + index.toString() as unknown as string);
    const svg = svgContainer.append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom);

    const g = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    const xScale = d3.scaleLinear()
        .domain([d3.min(dat, data => d3.min(data, d => d.x))!, d3.max(dat, data => d3.max(data, d => d.x))!])
        .range([0, width]);

    const yScale = d3.scaleLinear()
        .domain([d3.max(dat, data => d3.max(data, d => d.y))!, 0])
        .range([0, height]);

	const line = d3.line<SpectrumPoint>()
		.x(d => xScale(d.x))
		.y(d => yScale(d.y));

    let xAxis = g.append("g")
        .attr("stroke", "white")
        .attr("stroke-width", 1)
        .attr("transform", "translate(0," + height + ")");
    xAxis.call(d3.axisBottom(xScale).ticks(10));

    let yAxis = g.append("g")
        .attr("stroke", "white")
        .attr("stroke-width", 1);
    yAxis.call(d3.axisLeft(yScale).tickFormat(d3.format('.0e')));

	yAxis.selectAll("text")
		.style("font-size", "16px");

	xAxis.selectAll("text")
		.style("font-size", "16px");

    dat.forEach((data, i) => {
        const color = colors[i % colors.length];
        
        // Draw the sticks of the lollipops
        g.selectAll(`.stick-${i}`)
			.data(data)
			.enter()
			.append('line')
			.attr('class', `stick stick-${i}`)
            .attr('x1', d => xScale(d.x))
            .attr('y1', d => yScale(d.y))
            .attr('x2', d => xScale(d.x))
            .attr('y2', height)
            .attr('stroke', color)
            .attr('stroke-width', 3);
    });

	const lines: string[] = dat.map((_, i) => `.stick-${i}`);
	const checkboxesContainer = svg.append("g")
		.attr("transform", `translate(${width}, ${0})`);


	function updateVisibility(index: number, checked: boolean) {
		g.selectAll(lines[index]).attr('visibility', checked ? 'visible' : 'hidden');
	}

    function renderCheckboxes() {
		checkboxesContainer.selectAll('*').remove();

		dat.forEach((_, i) => {
			const color = colors[i % colors.length];
        	const checkboxGroup = checkboxesContainer.append("g")
				.attr("transform", `translate(${(-i * 160)-40}, ${10})`);

			const isChecked = g.select(`.stick-${i}`).attr('visibility') !== 'hidden';


			checkboxGroup.append("foreignObject")
				.attr("width", 15)
				.attr("height", 15)
				.attr("y", -5)
				.append("xhtml:input")
				.attr("type", "checkbox")
				.attr("checked", isChecked)
				.on("change", function() {
					const checked = d3.select(this).property("checked");
					updateVisibility(i, checked);
					if ((this! as HTMLInputElement).checked) {
						shouldIncludeArray[i] = true;
					} else {
						shouldIncludeArray[i] = false;
					}
					zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, SpectrumPoint>);
				});


			checkboxGroup.append("foreignObject")
				.attr("x", 60) // Setting this to a value so it's to the right of the checkbox and label.
				.attr("width", 70)  // Adjust width as per your requirements.
				.attr("height", 20) // Adjust height as per your requirements.
				.attr("y", -10)
				.append("xhtml:button")
				.attr("class", "ms2-compare-button")
				.text("Compare")
				.on("click", async function() {
					// Implement the compare functionality here.
					// You can use the current index "i" to determine which item to compare.
					
					let name = await invoke<string>("get_name_from_identifier_msms_tauri", {identifier:identifier})
					console.log(`Comparing item ${name} ${adduct} ${labels[i]}`);

					// Grab the div using its ID
					let myDiv: HTMLElement | null = document.getElementById("ms2-to-compare-compounds");

					if(myDiv) {
						let elementDiv: HTMLDivElement = document.createElement("div");
						elementDiv.className = "ms2-to-compare-button";
						elementDiv.id = `ms2-to-compare-button-${name}-${adduct}-${labels[i]}`

						let textSpan: HTMLSpanElement = document.createElement("span");
						textSpan.textContent = `${name} ${adduct} ${labels[i]}`;
						textSpan.className = "ms2-to-compare-textspan";
						textSpan.id = `${identifier}-${adduct}-${labels[i]}`;
						let iconSpan: HTMLSpanElement = document.createElement("span");
						iconSpan.className = "ms2-to-compare-iconspan";
						iconSpan.id = `ms2-to-compare-iconspan-${name}-${adduct}-${labels[i]}`;

						let icon: HTMLElement = document.createElement("ion-icon");
						icon.setAttribute("name", "close-outline");
						iconSpan.appendChild(icon);

						elementDiv.appendChild(textSpan);
						elementDiv.appendChild(iconSpan);

						myDiv.append(elementDiv);


						document.getElementById(`ms2-to-compare-iconspan-${name}-${adduct}-${labels[i]}`)!.addEventListener('click', function() {
							let targetDiv = document.getElementById(`ms2-to-compare-button-${name}-${adduct}-${labels[i]}`);
							if (targetDiv) {
								targetDiv.parentNode!.removeChild(targetDiv);
							}
						});

					}

				});

			checkboxGroup.append("text")
				.attr("x", 20)
				.attr("y", 10)
				.attr("fill", color)
				.text(`${labels[i]}`);
		});
	}

    function zoomed(event: d3.D3ZoomEvent<SVGSVGElement, SpectrumPoint>) {
        const transform = event.transform;
        const xNewScale = transform.rescaleX(xScale);
        const filteredData = dat.flatMap((data, index) => {
			if (shouldIncludeArray[index]) {
				return data.filter(d => d.x >= xNewScale.domain()[0] && d.x <= xNewScale.domain()[1]);
			} else {
				return []; // Exclude this data array
			}
		});
        const yMaxInView = d3.max(filteredData, d => d.y)!;
        const yNewScale = d3.scaleLinear().domain([yMaxInView, 0]).range([0, height]);

        xAxis.call(d3.axisBottom(xNewScale).ticks(10));
        yAxis.call(d3.axisLeft(yNewScale).tickFormat(d3.format('.0e')));

		yAxis.selectAll("text")
			.style("font-size", "12px");

		xAxis.selectAll("text")
			.style("font-size", "12px");
		
		g.selectAll<SVGLineElement, SpectrumPoint>('line.stick')
			.attr('x1', d => xNewScale(d.x))
			.attr('y1', d => yNewScale(d.y))
			.attr('x2', d => xNewScale(d.x))
			.attr('y2', height);
		
	}

    const zoom = d3.zoom<SVGSVGElement, unknown>()
        .scaleExtent([0.5, Infinity])
        .translateExtent([[0, 0], [width * 1.5, height]])
        .extent([[0, 0], [width, height]])
        .on('zoom', zoomed);

    svg.call(zoom);

    const resizeObserver = new ResizeObserver(() => {
		const newWidth = (svgContainer.node() as HTMLElement).getBoundingClientRect().width - margin.left - margin.right;
		const newHeight = (svgContainer.node() as HTMLElement).getBoundingClientRect().height - margin.top - margin.bottom;

		svg.attr('width', newWidth + margin.left + margin.right)
			.attr('height', newHeight + margin.top + margin.bottom);

		xScale.range([0, newWidth]);
		yScale.range([0, newHeight]);

		checkboxesContainer.attr("transform", `translate(${newWidth}, ${0})`);

		xAxis.attr('transform', 'translate(0,' + newHeight + ')')
			.call(d3.axisBottom(xScale).ticks(100));
		yAxis.call(d3.axisLeft(yScale).tickFormat(d3.format('.0e')));

		g.selectAll('.line').attr('d', d => line(d as SpectrumPoint[]));

		renderCheckboxes();
		zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, SpectrumPoint>);
	});

	resizeObserver.observe(svgContainer.node() as Element);
	zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, SpectrumPoint>);
}
