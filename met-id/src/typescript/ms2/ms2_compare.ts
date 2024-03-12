import * as d3 from 'd3';
import { invoke } from "@tauri-apps/api/tauri";
import { SpectrumPoint, SpectrumMap } from './ms2_draw_spectra';
import html2canvas from 'html2canvas';

let ExportListener: ((e: MouseEvent) => void) | null = null;
let ExportCsvListener: ((e: MouseEvent) => void) | null = null;

let ExportElement: HTMLButtonElement | null = null;
let ExportCsvElement: HTMLButtonElement | null = null;

function showPopup() {
    console.log("showing popup")
    const overlay = document.getElementById('ms2-overlay');
    const popup = document.getElementById('ms2-popup');

    overlay!.style.display = 'block';
    popup!.style.display = 'flex';
    
    const closePopupButton = document.getElementById("ms2-popup-cancel-button");
    closePopupButton?.addEventListener("click", () => hidePopup())

    const overlayDiv = document.getElementById("ms2-overlay");
    overlayDiv?.addEventListener("click", () => hidePopup());

    ExportElement = document.getElementById("ms2-popup-export-button") as HTMLButtonElement;
    ExportCsvElement = document.getElementById("ms2-popup-export-csv-button") as HTMLButtonElement;


    ExportListener = () => {
        save_msms_as_image()
    }
    ExportCsvListener = () => {
        // Select all span elements with the class 'ms2-to-compare-textspan'
        const spans = document.querySelectorAll('.ms2-to-compare-textspan');

        // Create an array to hold the text from each span
        const texts: string[] = [];
        const ids: string[] = [];
        const adducts: string[] = [];
        const cids: string[] = [];

        // Iterate over the NodeList of spans and add their text content to the array
        spans.forEach(span => {
            //console.log("the id is: ", span.id);
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

        drawMirrored(ids, adducts, cids, texts, true);
    }
    ExportElement?.addEventListener("click", ExportListener)
    ExportCsvElement?.addEventListener("click", ExportCsvListener)

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
    if (ExportListener) {
		ExportElement!.removeEventListener("click", ExportListener);
		ExportListener = null;
    }
    if (ExportCsvListener) {
		ExportCsvElement!.removeEventListener("click", ExportCsvListener);
		ExportCsvListener = null;
    }
}

function downloadCsv(csvContent, fileName) {
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = fileName;
    document.body.appendChild(anchor); // Required for Firefox
    anchor.click();
    document.body.removeChild(anchor);
    URL.revokeObjectURL(url);
}

async function drawMirrored(identifiers:string[], adducts:string[], cids:string[], texts:string[], exportCsv = false) {
    let ms2popupdiv = document.getElementById("ms2-popup-compare-names")!;
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
    const margin = { top: 20, right: 30, bottom: 50, left: 70 };
    
    const width = (svgContainer.node() as HTMLElement).getBoundingClientRect().width - margin.left - margin.right;
    const height = (svgContainer.node() as HTMLElement).getBoundingClientRect().height - margin.top - margin.bottom;

    const colors = ["red", "yellow", "blue", "purple", "orange", "green"];
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

        console.log(dat);
    // Handle regular and mirrored spectra separately
    for (let i = 0; i < dat.length; i++) {
        const data = dat[i];
        // Sort by height
        const sortedPeaks = data.sort((a, b) => b.y - a.y);

        const topPeaks: SpectrumPoint[] = [];
        sortedPeaks.forEach(peak => {
        // Check if the current peak is within 0.5 units on the x-axis of any already selected peak
        const isTooClose = topPeaks.some(selectedPeak => Math.abs(selectedPeak.x - peak.x) <= 0.5);
        // If it's not too close, and we haven't already selected 10 peaks, add it to the topPeaks
        if (!isTooClose && topPeaks.length < 10) {
            topPeaks.push(peak);
        }
        });

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
                        <div style="display: inline-block; margin-right: 30px; width: 12px; height: 10px; background-color: ${color}; border-radius: 10px;"></div>`;
        }

    
    ms2popupdiv!.innerHTML = HTMLInputs;
	//const _lines: string[] = dat.map((_, i) => `.stick-${i}`);
    
    function zoomed(event: d3.D3ZoomEvent<SVGSVGElement, unknown>) {
        const transform = event.transform;
        const xNewScale = transform.rescaleX(xScale);

        // Select and remove all circle elements
        d3.select('svg').selectAll('circle').remove();

        // Select and remove all text elements
        d3.select('svg').selectAll('text').remove();
        d3.select('svg').selectAll('g.tick').remove();
        
        //d3.select('svg').selectAll('g').remove();
        let yAxis = g.append("g")
            .attr("stroke", "white")
            .attr("stroke-width", 1);

        yAxis.selectAll("text")
            .style("font-size", "16px");
        
        
        
        for (let i = 0; i < dat.length; i++) {
            const mirroredMultiplier = i === 0 ? 1 : -1;
            const data = dat[i] as SpectrumPoint[];
            // 1. Sort by absolute y values, keeping the original sign
            const sortedByAbsoluteY = data.sort((a, b) => Math.abs(b.y) - Math.abs(a.y));

            // 2. Filter to ensure proximity rule on the x-axis
            const filteredForProximity = sortedByAbsoluteY.reduce((acc: SpectrumPoint[], curr) => {
            const tooClose = acc.some(item => Math.abs(item.x - curr.x) <= 0.5);
            if (!tooClose) acc.push(curr);
                return acc;
            }, []);

            // Limit to top 10 after filtering
            const top10Peaks: SpectrumPoint[] = filteredForProximity.slice(0, 10);
    
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

            top10Peaks.forEach(peak => {
                svg.append('text')
                    .attr('x', xNewScale(peak.x) + margin.left)
                    .attr('y', yNewScale(peak.y*mirroredMultiplier) - (12*mirroredMultiplier) + margin.top) // Adjust position above the circle
                    .attr('color', 'white')
                    .attr('text-anchor', 'middle')
                    .style('font-size', '12px')
                    .style('fill', 'white')
                    .text(`${peak.x.toFixed(2)}`);
                });
            }
        let xAxis = g.append("g")
            .attr("stroke", "white")
            .attr("stroke-width", 1)
            .attr("transform", "translate(0," + height + ")");
        xAxis.call(d3.axisBottom(xNewScale).ticks(10));

        xAxis.selectAll("text")
            .style("font-size", "16px");
    
    }
        

    const zoom = d3.zoom<SVGSVGElement, unknown>()
        .scaleExtent([0.5, Infinity])
        .translateExtent([[0, 0], [width * 1.3, height]])
        .extent([[0, 0], [width, height]])
        .on('zoom', zoomed);

    svg.call(zoom);

	zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, unknown>);

    function spectrumPointsToCsv() {
        // Calculate the max length of sub-arrays to define headers
        const maxLength = dat.reduce((max, arr) => Math.max(max, arr.length), 0);
        
        // Header
        let headers = '';
        for (let i = 0; i < dat.length; i++) {
            headers += `${texts[i]} m/z, ${texts[i]} Intensity,`;
        }
        headers = headers.slice(0, -1); // Remove the trailing comma
        let csvContent = `${headers}\n`;
    
        // Data rows
        for (let rowIndex = 0; rowIndex < maxLength; rowIndex++) {
            let row = '';
            for (let arrIndex = 0; arrIndex < dat.length; arrIndex++) {
                const point = dat[arrIndex][rowIndex];
                if (point) {
                    row += `${point.x}, ${point.y}, `;
                } else {
                    row += ',,'; // Empty placeholders if no point exists at this index
                }
            }
            row = row.slice(0, -1); // Remove the trailing comma
            csvContent += `${row}\n`;
        }
    
        // Assuming downloadCsv is defined elsewhere to handle the CSV download
        downloadCsv(csvContent, 'spectrumPoints.csv');
        zoomed({ transform: d3.zoomIdentity } as d3.D3ZoomEvent<SVGSVGElement, unknown>);
    }
    
    
    if (exportCsv) {
        spectrumPointsToCsv()
    }
}

function convertSvgToCanvas(svgElement, callback) {
    const svgData = new XMLSerializer().serializeToString(svgElement);
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
    const img = new Image();

    img.onload = () => {
        canvas.width = img.width;
        canvas.height = img.height;
        if (ctx) {
            ctx.drawImage(img, 0, 0);
            callback(canvas); // Pass the canvas to the callback
        }
    };

    const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
    const URL = window.URL || window.webkitURL || window;
    const blobURL = URL.createObjectURL(svgBlob);

    img.src = blobURL;
}

function saveCanvasAsImage(canvas: HTMLCanvasElement, filename: string): void {
    const imageURL = canvas.toDataURL('image/png').replace('image/png', 'image/octet-stream');
    const link = document.createElement('a');
    link.download = filename;
    link.href = imageURL;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

export function save_msms_as_image() {
    console.log("saving");

    let a = document.getElementById("ms2-popup-compare-names");
    // Step 1: Capture HTML content as canvas
    html2canvas(a as HTMLElement, {backgroundColor: '#2F2F2F'}).then(htmlCanvas => {
        // Assuming you have an SVG element for your D3 plot
        const svgElement = document.querySelector('svg') as SVGSVGElement;

        // Step 2: Convert SVG to canvas
        convertSvgToCanvas(svgElement, (svgCanvas) => {
            // Now you have both HTML content and SVG as canvas elements

            // Step 3: Combine both canvases
            const combinedCanvas = document.createElement('canvas');
            const ctx = combinedCanvas.getContext('2d');

            // Adjust combined canvas size to fit both
            combinedCanvas.width = Math.max(htmlCanvas.width, svgCanvas.width);
            combinedCanvas.height = htmlCanvas.height + svgCanvas.height;

            console.log(combinedCanvas)

            if (ctx) {
                ctx.fillStyle = '#2F2F2F'; // Set the background color
                ctx.fillRect(0, 0, combinedCanvas.width, combinedCanvas.height); // Fill the background
                // Draw HTML canvas content first
                ctx.drawImage(htmlCanvas, 50, 0);

                // Then draw SVG canvas below it
                ctx.drawImage(svgCanvas, 0, htmlCanvas.height);

                // Step 4: Save the combined canvas as an image
                saveCanvasAsImage(combinedCanvas, 'myPlotWithHtml.png');
            }
        });
    });
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
        //console.log("the id is: ", span.id);
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

