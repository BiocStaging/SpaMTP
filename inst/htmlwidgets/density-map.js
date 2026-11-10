// Development assistance: see inst/COPYRIGHTS.

function(el, x, json) {
    "use strict";
    const data = JSON.parse(json);
    const previous = document.getElementById(el.id + "-density");
    if (previous) previous.remove();
    const panel = document.createElement("section");
    panel.id = el.id + "-density";
    panel.style.cssText = "font:14px sans-serif;padding:16px;max-width:960px;margin:auto";
    el.parentNode.insertBefore(panel, el.nextSibling);
    const title = document.createElement("h2");
    title.textContent = "Spatial detection density";
    panel.appendChild(title);
    const controls = document.createElement("div");
    panel.appendChild(controls);
    const featureLabel = document.createElement("label");
    featureLabel.textContent = "Feature ";
    const select = document.createElement("select");
    data.features.forEach(function(feature, i) {
        const option = document.createElement("option");
        option.value = i;
        option.textContent = feature.feature + " (m/z " + feature.mz + ")";
        select.appendChild(option);
    });
    featureLabel.appendChild(select);
    controls.appendChild(featureLabel);
    const bandwidthLabel = document.createElement("label");
    bandwidthLabel.textContent = " Bandwidth (display pixels) ";
    const bandwidth = document.createElement("input");
    bandwidth.type = "number";
    bandwidth.min = "1";
    bandwidth.max = "100";
    bandwidth.value = "20";
    bandwidthLabel.appendChild(bandwidth);
    controls.appendChild(bandwidthLabel);
    const summary = document.createElement("p");
    panel.appendChild(summary);
    const svg = d3.select(panel).append("svg").attr("viewBox", "0 0 720 480")
        .attr("role", "img").attr("aria-label", "Spatial detections and kernel density contours");
    const scale = Math.min(660 / (data.xlim[1] - data.xlim[0] || 1),
        420 / (data.ylim[1] - data.ylim[0] || 1));
    const px = function(d) { return 30 + (d[0] - data.xlim[0]) * scale; };
    const py = function(d) { return 30 + (d[1] - data.ylim[0]) * scale; };
    function draw() {
        const i = Number(select.value);
        const feature = data.features[i];
        const points = data.locations[i];
        const width = Math.max(1, Math.min(100, Number(bandwidth.value) || 20));
        svg.selectAll("*").remove();
        summary.textContent = points.length + " detected pixels. " +
            (feature.annotation || "No annotation") + " " + (feature.identifiers || "");
        if (points.length) {
            const contours = d3.contourDensity().x(px).y(py).size([720, 480])
                .bandwidth(width).thresholds(15)(points);
            const maximum = d3.max(contours, function(d) { return d.value; }) || 1;
            const color = d3.scaleSequential(d3.interpolateViridis).domain([0, maximum]);
            svg.append("g").selectAll("path").data(contours).join("path")
                .attr("d", d3.geoPath()).attr("fill", function(d) { return color(d.value); });
            svg.append("g").selectAll("circle").data(points).join("circle")
                .attr("cx", px).attr("cy", py).attr("r", 1.5)
                .attr("fill", "black").attr("opacity", 0.35);
        }
    }
    select.addEventListener("change", draw);
    bandwidth.addEventListener("input", draw);
    el.on("plotly_click", function(event) {
        const point = event.points[0];
        if (Number.isInteger(point.customdata)) {
            select.value = String(point.customdata);
            draw();
        }
    });
    draw();
}
