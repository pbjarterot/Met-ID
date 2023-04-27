import * as d3 from 'd3';


export function draw() {
    console.log("Hello there");

    interface DataPoint {
    x: number;
    y: number;
    }

    const data: DataPoint[] = [
    { x: 0, y: 5 },
    { x: 1, y: 9 },
    { x: 2, y: 7 },
    { x: 3, y: 5 },
    { x: 4, y: 3 },
    { x: 5, y: 8 },
    ];

    const margin = { top: 20, right: 30, bottom: 30, left: 40 };
    const width = 600 - margin.left - margin.right;
    const height = 400 - margin.top - margin.bottom;

    const xScale = d3
    .scaleBand<number>()
    .range([0, width])
    .domain(data.map((d) => d.x))
    .padding(0.1);

    const yScale = d3
    .scaleLinear()
    .range([height, 0])
    .domain([0, d3.max(data, (d) => d.y)!]);

    const lineGenerator = d3
    .line<DataPoint>()
    .x((d) => xScale(d.x) ?? 0)
    .y((d) => yScale(d.y) ?? 0);

    const container = d3.select('#my-chart-container');
    console.log(container);
    const svg = container
        .append('svg')
        .attr('width', width + margin.left + margin.right)
        .attr('height', height + margin.top + margin.bottom)
        .append('g')
        .attr('transform', `translate(${margin.left},${margin.top})`);

    svg
        .append('path')
        .datum(data)
        .attr('fill', 'none')
        .attr('stroke', 'steelblue')
        .attr('stroke-width', 2)
        .attr('d', lineGenerator);
}
