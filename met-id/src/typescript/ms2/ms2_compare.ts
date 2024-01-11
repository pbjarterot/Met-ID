import * as d3 from 'd3';
import { invoke } from "@tauri-apps/api/tauri";
import { SpectrumPoint, SpectrumMap } from './ms2_draw_spectra';

function showPopup() {
    const overlay = document.getElementById('ms2-overlay');
    const popup = document.getElementById('ms2-popup');

    overlay!.style.display = 'block';
    popup!.style.display = 'flex';
    
    const closePopupButton = document.getElementById("ms2-popup-cancel-button");
    closePopupButton?.addEventListener("click", () => hidePopup())

    const overlayDiv = document.getElementById("ms2-overlay");
    overlayDiv?.addEventListener("click", () => hidePopup());


}

function hidePopup() {
    const overlay = document.getElementById('ms2-overlay');
    const popup = document.getElementById('ms2-popup');

    overlay!.style.display = 'none';
    popup!.style.display = 'none';

    const plot = document.getElementById('ms2-popup-plot-svg') as HTMLDivElement;
    if (plot) {
        plot.parentNode!.removeChild(plot);
    }

}

async function drawMirrored(identifiers:string[], adducts:string[], cids:string[], texts:string[]) {
    let ms2popupdiv = document.getElementById("ms2-popup-compare-names")!
    let HTMLInputs = "";

    


    const fetchData = async (identifiers: string[], adducts: string[], cids:string[]): Promise<{ keys: string[], values: SpectrumPoint[][] }> => {
        const keys: string[] = [];
        const values: SpectrumPoint[][] = [];
        for (let i = 0; i < identifiers.length; i++) {
            const result = await invoke<string>('get_msms_spectra_tauri', { identifier:identifiers[i], adduct:adducts[i], cid:cids[i] });
            const resultMap = JSON.parse(result) as SpectrumMap;
            for (const key in resultMap) {
                keys.push(key);
                values.push(resultMap[key]);
            }
            
        };
        return { keys, values };
    };
    const svgContainer = d3.select("#ms2-popup-plot" as unknown as string);
    const margin = { top: 80, right: 30, bottom: 30, left: 70 };
    
    const width = (svgContainer.node() as HTMLElement).getBoundingClientRect().width - margin.left - margin.right;
    const height = (svgContainer.node() as HTMLElement).getBoundingClientRect().height - margin.top - margin.bottom;

    const colors = ["red", "yellow", "blue", "purple", "orange", "yellow"];
    const {keys: _labels, values: dat} = await fetchData(identifiers, adducts, cids);
    
    const svg = svgContainer.append("svg")
        .attr("id", "ms2-popup-plot-svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom);

    const g = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let yMax = d3.max(dat, data => d3.max(data, d => d.y));
    yMax = yMax || 0;  // Ensure yMax is a number
    // Adjust yScale for mirrored spectra
    const yScale = d3.scaleLinear()
        .domain([-yMax, 0, yMax])
        .range([0, height/2, height]);

    const xScale = d3.scaleLinear()
        .domain([d3.min(dat, data => d3.min(data, d => d.x))!, d3.max(dat, data => d3.max(data, d => d.x))!])
        .range([0, width]);
    

    let yAxis = g.append("g")
        .attr("stroke", "white")
        .attr("stroke-width", 1);

	yAxis.selectAll("text")
		.style("font-size", "16px");
    let xAxis = g.append("g")
        .attr("stroke", "white")
        .attr("stroke-width", 1)
        .attr("transform", "translate(0," + height + ")");
    xAxis.call(d3.axisBottom(xScale).ticks(10));

	xAxis.selectAll("text")
		.style("font-size", "16px");

    // Handle regular and mirrored spectra separately
    for (let i = 0; i < dat.length; i++) {
        const data = dat[i];
        const color = colors[i % colors.length];
        const mirroredMultiplier = i === 0 ? 1 : -1;  // Regular for first array, mirrored for second
        // Draw the sticks of the lollipops
        g.selectAll(`.stick-${i}`)
            .data(data)
            .enter()
            .append('line')
            .attr('class', `stick stick-${i}`)
            .attr('x1', d => xScale(d.x))
            .attr('y1', yScale(0)) // Start from the center (0) for mirrored spectra
            .attr('x2', d => xScale(d.x))
            .attr('y2', d => yScale(d.y* mirroredMultiplier)) // Apply mirrored effect
            .attr('stroke', color)
            .attr('stroke-width', 3);

            HTMLInputs += ` 
                            <p style="display: inline-block;">${texts[i]}</p>
                            <div style="display: inline-block; margin-right: 30px; width: 12px; height: 10px; background-color: ${color}; border-radius: 10px;"></div>`
    }

    
    ms2popupdiv!.innerHTML = HTMLInputs;
	//const _lines: string[] = dat.map((_, i) => `.stick-${i}`);
    
    function zoomed(event: d3.D3ZoomEvent<SVGSVGElement, unknown>) {
        const transform = event.transform;
        const xNewScale = transform.rescaleX(xScale);
    
        xAxis.call(d3.axisBottom(xNewScale).ticks(10));
        xAxis.selectAll("text")
		.style("font-size", "16px");
        yAxis.selectAll("text")
		.style("font-size", "16px");
    
        for (let i = 0; i < dat.length; i++) {
            const mirroredMultiplier = i === 0 ? 1 : -1;
            const data = dat[i] as SpectrumPoint[];
    
            const filteredData = data.filter(d => d.x >= xNewScale.domain()[0] && d.x <= xNewScale.domain()[1]);
            const yMaxInView = d3.max(filteredData, d => Math.abs(d.y)*1.2)!;
    
            const yNewScale = d3.scaleLinear()
                                .domain([-yMaxInView, yMaxInView])
                                .range([height, 0]);
    
            yAxis.call(d3.axisLeft(yNewScale).tickFormat(d3.format('.0e')));
    
            g.selectAll(`line.stick-${i}`)
                .attr('x1', d => xNewScale((d as SpectrumPoint).x))
                .attr('y1', yScale(0))
                .attr('x2', d => xNewScale((d as SpectrumPoint).x))
                .attr('y2', d => yNewScale((d as SpectrumPoint).y * mirroredMultiplier));
        }
    }

    const zoom = d3.zoom<SVGSVGElement, unknown>()
        .scaleExtent([0.5, Infinity])
        .translateExtent([[0, 0], [width * 1.5, height]])
        .extent([[0, 0], [width, height]])
        .on('zoom', zoomed);

    svg.call(zoom);

	zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, unknown>);
    
}



export function compare_msms() {
    showPopup();
    // Select all span elements with the class 'ms2-to-compare-textspan'
    const spans = document.querySelectorAll('.ms2-to-compare-textspan');

    // Create an array to hold the text from each span
    const texts: string[] = [];
    const ids: string[] = [];
    const adducts: string[] = [];
    const cids: string[] = [];

    // Iterate over the NodeList of spans and add their text content to the array
    spans.forEach(span => {
        console.log("the id is: ", span.id);
        if (span.textContent) {
            if (span.textContent !== "Compare!") {
                texts.push(span.textContent);
                
                let words = span.id.split("-"); // Splits the string into an array of words
                ids.push(words[0]);
                adducts.push(words[1]); // Selects the last word
                cids.push(words[2]); // Selects the last word
            }
            
        }
        
    });

    drawMirrored(ids, adducts, cids, texts);

}