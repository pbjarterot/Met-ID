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



async function drawMirrored(_index: number, identifier:string, adduct:string) {
    const fetchData = async (identifier: string, adduct: string): Promise<{ keys: string[], values: SpectrumPoint[][] }> => {
        const result = await invoke<string>('get_msms_spectra', { identifier, adduct });
        const resultMap = JSON.parse(result) as SpectrumMap;

        const keys: string[] = [];
        const values: SpectrumPoint[][] = [];

        for (const key in resultMap) {
            keys.push(key);
            values.push(resultMap[key]);
        }

        return { keys, values };
    };
    const svgContainer = d3.select("#ms2-popup-plot" as unknown as string);
    const margin = { top: 100, right: 30, bottom: 30, left: 70 };
    
    const width = (svgContainer.node() as HTMLElement).getBoundingClientRect().width - margin.left - margin.right;
    const height = (svgContainer.node() as HTMLElement).getBoundingClientRect().height - margin.top - margin.bottom;

    const colors = ["red", "green", "blue", "purple", "orange", "yellow"];
    const {keys: _labels, values: dat} = await fetchData(identifier, adduct);
    
    const svg = svgContainer.append("svg")
        .attr("id", "ms2-popup-plot-svg")
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
    /*
	const _line = d3.line<SpectrumPoint>()
		.x(d => xScale(d.x))
		.y(d => yScale(d.y));
    */

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

	//const _lines: string[] = dat.map((_, i) => `.stick-${i}`);




    function zoomed(event: d3.D3ZoomEvent<SVGSVGElement, SpectrumPoint>) {
        const transform = event.transform;
        const xNewScale = transform.rescaleX(xScale);
        const filteredData = dat.flatMap(data => data.filter(d => d.x >= xNewScale.domain()[0] && d.x <= xNewScale.domain()[1]));
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

	zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, SpectrumPoint>);
}



export function compare_msms() {
    showPopup();

    drawMirrored(0, "HMDB0000763", "1A");



}