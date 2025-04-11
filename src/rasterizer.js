"use strict";

const sizes = {
    width: window.innerWidth,
    height: window.innerHeight
}
const BLUE = {r: 0, g: 0, b: 255};
const RED = {r: 255, g: 0, b: 0};
const GREEN = {r: 0, g: 255, b: 0};

const d = 1;

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

class Vertex {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

class Triangle {
    constructor(v, color) {
        this.v = v;
        this.color = color;
    }
}

class Model {
    constructor(name, vertices, triangles) {
        this.name = name;
        this.vertices = vertices;
        this.triangles = triangles;
    }
}

class Transform {
    constructor(scale, rotation, translation) {
        this.scale = scale;
        this.rotation = rotation;
        this.translation = translation;
    }
}


class Instance {
    constructor(model, transform) {
        this.model = model;
        this.transform = transform;
    }
}

class Camera {
    constructor(position, rotation) {
        this.position = position;
        this.rotation = rotation;
    }
}

class Vec4 {
    constructor(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
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

const viewportToCanvas = (x, y) => {
    return {x: x * sizes.width / 1,y: y * sizes.height / 1}
}

const projectVertex = (v) => {
    return viewportToCanvas(v.x * d / v.z, v.y * d / v.z)
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

const drawWireframeTriangle = (p0, p1, p2, color) => {
    drawLine(p0, p1, color);
    drawLine(p1, p2, color);
    drawLine(p2, p0, color);
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

const renderTriangle = (triangle, projected) => {
    drawWireframeTriangle(projected[triangle.v[0]], projected[triangle.v[1]], projected[triangle.v[2]], triangle.color);
}

function makeOYRotationMatrix(rotation) {
    const degree = rotation.y;
    const cos = Math.cos(degree*Math.PI/180);
    const sin = Math.sin(degree*Math.PI/180);
    return [
        [cos, 0, -sin, 0],
        [0, 1, 0, 0],
        [sin, 0, cos, 0],
        [0, 0, 0, 1]
    ];
}

function makeTranslationMatrix(translation) {
    return [
        [1, 0, 0, translation.x],
        [0, 1, 0, translation.y],
        [0, 0, 1, translation.z],
        [0, 0, 0, 1]
    ];
}

function makeScaleMatrix(scale) {
    return [
        [scale.x, 0, 0, 0],
        [0, scale.y, 0, 0],
        [0, 0, scale.z, 0],
        [0, 0, 0, 1]
    ];
}

function multiplyMatrices(a, b) {
    const result = [];
    for (let i = 0; i < 4; i++) {
        result[i] = [];
        for (let j = 0; j < 4; j++) {
            let sum = 0;
            for (let k = 0; k < 4; k++) {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

// multiply matrix and 4d vector
function multiplyMatrixVector(matrix, vector) {
    const result = [];
    let vec = [vector.x, vector.y, vector.z, vector.w];

    for (let i = 0; i < 4; i++) {
        let sum = 0;
        for (let j = 0; j < 4; j++) {
            sum += matrix[i][j] * vec[j];
        }
        result[i] = sum;
    }
    return new Vec4(result[0], result[1], result[2], result[3]);
}

function transposeMatrix(matrix) {
    const result = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]];
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

const applyTransformation = (vertex, transform) => {
    const scaled = {
        x: vertex.x * transform.scale.x,
        y: vertex.y * transform.scale.y,
        z: vertex.z * transform.scale.z
    }
    
    const rotated = rotate(scaled, transform.rotation);
    const translated = translate(rotrated, transform.translation);
    return translated;
    
}

const renderObject = (instance, transformMatrix) => {
    let projected = []
    // console.log(instance[0]);
    const instancedModel = instance.model;
    const transformation = instance.transform;
    for (const vertex of instancedModel.vertices) {
        // translate vertex by {x: -1.5, y:0, z:7}
        // vertex.x += position.x;
        // vertex.y += position.y;
        // vertex.z += position.z;
        let projectedVertex = projectVertex(multiplyMatrixVector(transformMatrix, new Vec4(vertex.x, vertex.y, vertex.z, 1)))
        // console.log(multiplyMatrixVector(vertex))
        console.log(projectedVertex);
        projected.push(projectedVertex)
    }
    // console.log(projected)

    for (const triangle of instancedModel.triangles) {
        // console.log(triangle.v1, triangle.v2, triangle.v3)
        renderTriangle(triangle, projected)
    }
}

const makeCameraMatrix = (position, rotation) => {
    const translationMatrix = makeTranslationMatrix(position);
    const rotationMatrix = makeOYRotationMatrix(rotation);
    const cameraMatrix = multiplyMatrices(translationMatrix, rotationMatrix);
    return cameraMatrix;
}

const renderScene = () => {
    const cameraMatrix = makeCameraMatrix(camera.position, camera.rotation);
    // console.log(cameraMatrix);
    for (const cube of cubeInstance) {
        const transformMatrix = multiplyMatrices(
            makeTranslationMatrix(cube.transform.translation), 
            multiplyMatrices(
                makeScaleMatrix(cube.transform.scale), 
                makeOYRotationMatrix(cube.transform.rotation)
            )
        );
        
        
        const matrix = multiplyMatrices(cameraMatrix, transformMatrix);
        renderObject(cube, matrix);
    }
}

let camera = new Camera(new Vertex(-3, 1,2), {x: 0, y: -30, z: 0})

const v0 = new Vertex(1, 1, 1)
const v1 = new Vertex(-1, 1, 1)
const v2 = new Vertex(-1, -1, 1)
const v3 = new Vertex(1, -1, 1)

const v4 = new Vertex(1, 1, -1)
const v5 = new Vertex(-1, 1, -1)
const v6 = new Vertex(-1, -1, -1)
const v7 = new Vertex(1, -1, -1)

const vertices = [v0, v1, v2, v3, v4, v5, v6, v7]

// Triangle index for cube
const triangle_index = [
    new Triangle([0,1,2], RED),
    new Triangle([0,2,3], RED),
    new Triangle([4,0,3], GREEN),
    new Triangle([4,3,7], GREEN),
    new Triangle([5,4,7], BLUE),
    new Triangle([5,7,6], BLUE),
    new Triangle([1,5,6], GREEN),
    new Triangle([1,6,2], GREEN),
    new Triangle([4,5,1], BLUE),
    new Triangle([4,1,0], BLUE),
    new Triangle([2,6,7], RED),
    new Triangle([2,7,3], RED)
]

const cubeModel = new Model("cube", vertices, triangle_index)
const cubeA = new Instance(cubeModel, new Transform({x: 1.2,y: 1.2,z: 1.2}, {x: 0, y: Math.PI/4, z: 0}, {x: -1.5, y: 0, z: 7}))
const cubeB = new Instance(cubeModel, new Transform({x: 1.2,y: 1.2,z: 1.2}, {x: 0, y: Math.PI/4, z: 0}, {x: 1.25, y: 2, z: 7.5}))
const cubeInstance = [cubeA, cubeB]

renderScene(cubeInstance)

// create cube
// drawLine(projectVertex(vAf), projectVertex(vBf), BLUE);
// drawLine(projectVertex(vBf), projectVertex(vCf), BLUE);
// drawLine(projectVertex(vCf), projectVertex(vDf), BLUE);
// drawLine(projectVertex(vDf), projectVertex(vAf), BLUE);

// drawLine(projectVertex(vAb), projectVertex(vBb), RED);
// drawLine(projectVertex(vBb), projectVertex(vCb), RED);
// drawLine(projectVertex(vCb), projectVertex(vDb), RED);
// drawLine(projectVertex(vDb), projectVertex(vAb), RED);

// drawLine(projectVertex(vAf), projectVertex(vAb), GREEN);
// drawLine(projectVertex(vBf), projectVertex(vBb), GREEN);
// drawLine(projectVertex(vCf), projectVertex(vCb), GREEN);
// drawLine(projectVertex(vDf), projectVertex(vDb), GREEN);


// Draw your lines or other primitives here
// drawLine(P0, P1, {r: 255, g: 0, b: 0});
// drawLine(P1, P2, {r: 255, g: 0, b: 0});
// drawLine(P2, P0, {r: 255, g: 0, b: 0});
// drawShadedTriangle(P0, P1, P2, {r:0, g:255, b:0})
updateCanvas();