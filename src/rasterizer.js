"use strict";

const sizes = {
    width: window.innerWidth,
    height: window.innerHeight
}


const canvas = document.getElementById("canvas");
if (!canvas) {
    throw new Error("Canvas element not found");
}

const ctx = canvas.getContext("2d");
if (!ctx) {
    throw new Error("2D context not supported");
}

let canvas_buffer = ctx.getImageData(0, 0, sizes.width, sizes.height);

canvas.width = sizes.width;
canvas.height = sizes.height;

const putPixel = (x, y, color) => {
    x = canvas.width/2 + (x | 0);
    y = canvas.height/2 - (y | 0) - 1;

    if (x < 0 || x >= sizes.width || y < 0 || y >= sizes.height) {
        return;
    }

    // Calculate offset in the buffer
    const offset = 4 * (x + y * sizes.width);
    
    // Set pixel data with clamped color values
    canvas_buffer.data[offset] = color.r;
    canvas_buffer.data[offset + 1] = color.g;
    canvas_buffer.data[offset + 2] = color.b;
    canvas_buffer.data[offset + 3] = 255; // Alpha channel

}

const updateCanvas = () => {
    ctx.putImageData(canvas_buffer, 0, 0);
}

class Point {
    constructor(x, y, h) {
        this.x = x;
        this.y = y;
        this.h = h
    }
}


const P0 = new Point(-200, -250, 0.3);
const P1 = new Point(200, 50, 0.1);
const P2 = new Point(20, 250, 1.0);

const swapa = (a, b) => {
    let temp = a;
    a = b;
    b = temp;
    return [a, b];
}



const interpolate = (i0, d0, i1, d1) => {
    if (i0 == i1) {
        return [d0];
    }
    let values = []
    const a = (d1 - d0) / (i1 - i0);
    let d = d0;
    for (let i = i0; i <= i1; i++) {
        values.push(d);
        d += a
    }
    return values;
}

const drawLine = (p0, p1, color) => {
    let dx = p1.x - p0.x 
    let dy = p1.y - p0.y;
    if (Math.abs(dx) > Math.abs(dy)) {
        if (dx < 0) {
            // swap p0 and p1
            [p0, p1] = swapa(p0, p1);
        }

        let ys = interpolate(p0.x, p0.y, p1.x, p1.y);
        for (let x = p0.x; x <= p1.x; x++) {
            putPixel(x, ys[x - p0.x | 0], color);
        }
    } else {
        if (dy < 0) {
            [p0, p1] = swapa(p0, p1);
        }

        let xs = interpolate(p0.y, p0.x, p1.y, p1.x);
        for (let y = p0.y; y <= p1.y; y++) {
            putPixel(xs[y - p0.y | 0], y, color);
        }
    }
}

const drawFilledTriangle = (p0, p1, p2, color) => {
    // sort the points so that y0<=y1<=y2
    if (p1.y < p0.y) {
        [p0, p1] = swapa(p0, p1);
    }
    if (p2.y < p0.y) {
        [p0, p2] = swapa(p1, p2);
    }
    if (p2.y < p1.y) {
        [p1, p2] = swapa(p0, p1);
    }

    // compute the x coordinates of triangle edges
    let x01 = interpolate(p0.y, p0.x, p1.y, p1.x);
    let x12 = interpolate(p1.y, p1.x, p2.y, p2.x);
    let x02 = interpolate(p0.y, p0.x, p2.y, p2.x);

    // remove last eleement of x01 and concatenate the x01 and x12 arrays
    x01.pop();
    const x012 = x01.concat(x12);

    // Determining left and right sides
    let left = [];
    let right = [];
    const m = Math.floor(x012.length / 2);
    if (x02[m] < x012[m]) {
        left = x02;
        right = x012;
    } else {
        left = x012;
        right = x02;
    }

    // Draw horizontal line
    for (let y = p0.y; y <= p2.y; y++) {
        for (let x = left[y - p0.y]; x <= right[y - p0.y]; x++) {
            putPixel(x, y, color);
        }
    }
}

let h0 = 0.3
let h1 = 0.1
let h2 = 1.0

const drawShadedTriangle = (p0, p1, p2, color) => {
    // sort the points so that y0<=y1<=y2
    if (p1.y < p0.y) {
        [p0, p1] = swapa(p0, p1);
    }
    if (p2.y < p0.y) {
        [p0, p2] = swapa(p1, p2);
    }
    if (p2.y < p1.y) {
        [p1, p2] = swapa(p0, p1);
    }
    
    let x01 = interpolate(p0.y, p0.x, p1.y, p1.x);
    let h01 = interpolate(p0.y, p0.h, p1.y, p1.h);
    
    let x12 = interpolate(p1.y, p1.x, p2.y, p2.x);
    let h12 = interpolate(p1.y, p1.h, p2.y, p2.h);
    
    let x02 = interpolate(p0.y, p0.x, p2.y, p2.x);
    let h02 = interpolate(p0.y, p0.h, p2.y, p2.h);

    // remove last element
    x01.pop();
    let x012 = x01.concat(x12);

    h01.pop();
    let h012 = h01.concat(h12);

    // Determining left and right sides
    let x_left = [];
    let h_left = [];
    let x_right = [];
    let h_right = [];

    const m = Math.floor(x012.length / 2) | 0;
    if (x02[m] < x012[m]) {
        x_left = x02;
        h_left = h02;

        x_right = x012;
        h_right = h012;
    } else {
        x_left = x012;
        h_left = h012;

        x_right = x02;
        h_right = h02;
    }

    // Draw horizontal segments
    for (let y = p0.y; y <= p2.y; y++) {
        let x_l = x_left[y - p0.y] | 0
        let x_r = x_right[y - p0.y] | 0
        
        const h_segment = interpolate(x_l, h_left[y - p0.y], x_r, h_right[y - p0.y]);
         
        for (let x = x_l; x <= x_r; x++) {
            const shaded_color = {
                r: color.r * h_segment[x - x_l],
                g: color.g * h_segment[x - x_l],
                b: color.b * h_segment[x - x_l]
            }
            putPixel(
                x, 
                y, 
                shaded_color
            );
        }
    }
}

// Draw your lines or other primitives here
drawLine(P0, P1, {r: 255, g: 0, b: 0});
drawLine(P1, P2, {r: 255, g: 0, b: 0});
drawLine(P2, P0, {r: 255, g: 0, b: 0});
drawShadedTriangle(P0, P1, P2, {r:0, g:255, b:0})
updateCanvas();