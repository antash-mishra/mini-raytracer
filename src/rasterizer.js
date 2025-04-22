"use strict";

// size for the canvas
const sizes = {
    width: window.innerWidth,
    height: window.innerHeight
}
const RED = {r: 255, g:0, b:0};
const GREEN = {r:0, g:255, b:0};
const BLUE = {r:0, g:0, b:255};
const YELLOW = {r:255, g:255, b:0};
const PURPLE = {r:255, g:0, b:255};
const CYAN = {r:0, g:255, b:255};

// distnace between viewport and camera in z axis 
const d = 1;

// canvas
const canvas = document.getElementById("canvas");
canvas.width = sizes.width
canvas.height = sizes.height

if (!canvas) {
    throw new Error("Canvas element not found");
}


let depth_buffer = Array();
depth_buffer.length = canvas.width * canvas.height;

function UpdateDepthBufferIfCloser(x, y, inv_z) {
    x = canvas.width/2 + (x | 0);
    y = canvas.height/2 - (y | 0) - 1;
  
    if (x < 0 || x >= canvas.width || y < 0 || y >= canvas.height) {
      return false;
    }
    
    let offset = x + canvas.width*y;
  if (depth_buffer[offset] == undefined || depth_buffer[offset] < inv_z) {
    depth_buffer[offset] = inv_z;
    return true;
  }
  return false;
}


const ctx = canvas.getContext("2d");
if (!ctx) {
    throw new Error("2D context not supported");
}

let canvas_buffer = ctx.getImageData(0, 0, sizes.width, sizes.height);


// set cavas color on the x and y axis
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

// class to represent a point in 2D space where h is the color intensity at that pixel
class Point {
    constructor(x, y, h) {
        this.x = x;
        this.y = y;
        this.h = h
    }
}

// class to represent a vertex in 3D space
class Vertex {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

// class to represent a triangle in 3D space where v is the index of the vertices in the model
// and color is the color of the triangle
class Triangle {
    constructor(v, color) {
        this.v = v;
        this.color = color;
    }
}

class Sphere {
    constructor(divs, color) {
        this.divs = divs;
        this.color = color;
    }
}

// class to represent a model in 3D space where vertices is the array of vertices in the model
// triangles is the array of triangles in the model, bounding_center is the center of the bounding sphere
// and bounding_radius is the radius of the bounding sphere
class Model {
    constructor(name, vertices, triangles, bounding_center, bounding_radius) {
        this.name = name;
        this.vertices = vertices;
        this.triangles = triangles;
        this.bounding_center = bounding_center;
        this.bounding_radius = bounding_radius;
    }
}

// class to represent a transform in 3D space where scale is the scale of the model, rotation is the rotation of the model
// and translation is the translation of the model
class Transform {
    constructor(scale, rotation, translation) {
        this.scale = scale;
        this.rotation = rotation;
        this.translation = translation;
    }
}

// class to represent an instance of a model in 3D space where model is the model and transform is the transform of the model
class Instance {
    constructor(model, transform) {
        this.model = model;
        this.transform = transform;
    }
}

// class to represent camera in 3D space where position is the position of the camera and rotation is the rotation of the camera
class Camera {
    constructor(position, rotation) {
        this.position = position;
        this.rotation = rotation;
    }
}

// class to represent a 4D vector where x, y, z are the coordinates of the vector
// and w is the homogeneous coordinate that represent vector being 1 and point being 0
class Vec4 {
    constructor(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }
}

// class to represent a plane in 3D space where normal is the normal vector of the plane
// and d is the distance from the origin to the plane
class Plane {
    constructor(normal, d) {
        this.normal = normal;
        this.d = d;
    }
}

class Light {
    constructor(type, intensity, vector) {
        this.type = type;
        this.intensity = intensity;
        this.vector = vector;
    }
}

const P0 = new Point(-200, -250, 0.3);
const P1 = new Point(200, 50, 0.1);
const P2 = new Point(20, 250, 1.0);

// Function to swap two variables
const swapa = (a, b) => {
    let temp = a;
    a = b;
    b = temp;
    return [a, b];
}

// Function to convert viewport coordinates to canvas coordinates
const viewportToCanvas = (x, y) => {
    return { x: x * sizes.width / 1 | 0,y: y * sizes.height / 1 | 0 }
}

// Function to project a 3D vertex to 2D canvas coordinates
const projectVertex = (v) => {
    return viewportToCanvas(v.x * d / v.z, v.y * d / v.z)
}

function GenerateSphere(divs, color) {
    let vertices = [];
    let triangles = [];

    let delta_angle = 2.0 * Math.PI / divs;
    // Generate vertices and normals
    for (let d=0; d<divs+1; d++) {
        let y = (2.0 / divs) * (d - divs / 2);
        let radius = Math.sqrt(1.0 - y * y);

        for (let i=0; i<divs; i++) {
            let vertex = new Vertex(radius * Math.cos(delta_angle * i), y, radius * Math.sin(delta_angle * i));
            vertices.push(vertex);
        }
    }

    // Generate triangles
    for (let d=0; d<divs; d++) {
        for (let i=0; i<divs; i++) {
            let i0 = d * divs + i;
            let i1 = (d+1) * divs + (i+1) % divs;
            let i2 = divs*d + (i+1)%divs;
            let tri0 = [i0, i1, i2];
            let tri1 = [i0, i0+divs, i1];
            triangles.push(new Triangle(tri0, color));
            triangles.push(new Triangle(tri1, color));
        }
    
    }
    console.log(vertices, triangles)
    return new Model("sphere", vertices, triangles, new Vertex(0, 0, 0), 1);
}

const dot = (a, b) => {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

const magnitude = (v) => {
    const result = Math.sqrt(dot(v, v))
    return result;
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
        [p0, p2] = swapa(p0, p2);
    }
    if (p2.y < p1.y) {
        [p1, p2] = swapa(p1, p2);
    }

    console.log(p0, p1, p2)

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


const sortVertices = (vertex_indexes,  projected) => {
    let indexes = [0, 1, 2]

    if (projected[vertex_indexes[indexes[1]]].y < projected[vertex_indexes[indexes[0]]].y) {
        let swap = indexes[0]; indexes[0] = indexes[1]; indexes[1] = swap; 
    }
    if (projected[vertex_indexes[indexes[2]]].y < projected[vertex_indexes[indexes[0]]].y) {
        let swap = indexes[0]; indexes[0] = indexes[2]; indexes[2] = swap; 
    }
    if (projected[vertex_indexes[indexes[2]]].y < projected[vertex_indexes[indexes[1]]].y) {
        let swap = indexes[1]; indexes[1] = indexes[2]; indexes[2] = swap; 
    }
    return indexes;
}

const computeIllumination = (vertex, normal, camera, lights) => {
    let illumination = 0;

    for (const light of lights) {
        if (light.type == LT_AMBIENT) {
            illumination += light.intensity;
            continue;
        }

        let vl;
        if (light.type == LT_DIRECTIONAL) {
            let cameraMatrix = transposeMatrix(makeOYRotationMatrix(camera.rotation));

            let rotatedLight = multiplyMatrixVector(cameraMatrix, new Vec4(light.vector.x, light.vector.y, light.vector.z, 1));
            vl = rotatedLight;
        }
        else if (light.type == LT_POINT) {
            let cameraMatrix = multiplyMatrices(
                makeOYRotationMatrix(camera.rotation), 
                makeTranslationMatrix(camera.position)
            );
            let transformed_light = multiplyMatrixVector(
                cameraMatrix, 
                new Vec4(light.vector.x, light.vector.y, light.vector.z, 1)
            );
            vl = new Vec4(
                transformed_light.x - vertex.x,
                transformed_light.y - vertex.y,
                transformed_light.z - vertex.z,
                1
            )
        }
        
        const n_dot_l = dot(normal, vl);

        if (LM_DIFFUSE) {
            if (n_dot_l > 0) {
                illumination += light.intensity * n_dot_l / (magnitude(normal) * magnitude(vl));
            }
        }

        if (LM_SPECULAR) {
            let reflected = { x: 2 * n_dot_l * normal.x - vl.x, y: 2 * n_dot_l * normal.y - vl.y, z: 2 * n_dot_l * normal.z - vl.z };
            let view = { x: camera.position.x - vertex.x, y: camera.position.y - vertex.y, z: camera.position.z - vertex.z };
            const cos_beta = dot(reflected, view) / (reflected.length * view.length);
            if (cos_beta > 0) {
                let specular = 50
                illumination += light.intensity * Math.pow(cos_beta, light.specular);
            }
        }
    } 
    return illumination;
}

const edgeInterpolate = (y0, v0, y1, v1, y2, v2) => {
    let v01 = interpolate(y0, v0, y1, v1);
    let v12 = interpolate(y1, v1, y2, v2);
    let v02 = interpolate(y0, v0, y2, v2);
    v01.pop();
    let v012 = v01.concat(v12);
    return [v02, v012];
}

const renderTriangle = (triangle, projected, vertices) => {
    // sort the vertices from the projected vertices
    let [i0, i1, i2] = sortVertices(triangle.v, projected);

    // Vertcies of triangle
    let v0 = vertices[triangle.v[i0]];
    let v1 = vertices[triangle.v[i1]];
    let v2 = vertices[triangle.v[i2]];


    // create normal for this triangle using the vertices
    // Using unsorted vertices is ikportant as it will help in finding back and front face properly
    // first vector
    let a = new Vertex(
        vertices[triangle.v[1]].x - vertices[triangle.v[0]].x, 
        vertices[triangle.v[1]].y - vertices[triangle.v[0]].y, 
        vertices[triangle.v[1]].z - vertices[triangle.v[0]].z
    );

    // second vector
    let b = new Vertex(
        vertices[triangle.v[2]].x - vertices[triangle.v[0]].x, 
        vertices[triangle.v[2]].y - vertices[triangle.v[0]].y, 
        vertices[triangle.v[2]].z - vertices[triangle.v[0]].z);

    // cross product
    let normal = new Vertex(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
    // normalize the normal vector
    let length = Math.sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    normal.x /= length;
    normal.y /= length;
    normal.z /= length;



    // Checking if angle between the normal and the camera is less than 90 degrees or dot product is greater than 0
    const vertex_to_camera = vertices[triangle.v[0]] - camera.position;
    const dot_product = normal.x * vertex_to_camera.x + normal.y * vertex_to_camera.y + normal.z * vertex_to_camera.z;
    if (dot_product < 0) {
        // console.log("backface culling")
        return;
    }



    let p0 = projected[triangle.v[i0]];
    let p1 = projected[triangle.v[i1]];
    let p2 = projected[triangle.v[i2]];
  
    // calculate attribute values at the edges
    let [x02, x012] = edgeInterpolate(p0.y, p0.x, p1.y, p1.x, p2.y, p2.x); // for drawing color
    let [z02, z012] = edgeInterpolate(p0.y, 1.0/v0.z, p1.y, 1.0/v1.z, p2.y, 1.0/v2.z); // for depth buffer calculation


    // Compute Flat Shading
    // calculate center of triangle
    let center = new Vertex(
        (v0.x + v1.x + v2.x) / 3,
        (v0.y + v1.y + v2.y) / 3,
        (v0.z + v1.z + v2.z) / 3
    );
    let intensity = computeIllumination(center, normal, camera, lights);

    
    // Determining left and right sides
    let m = (x02.length / 2) | 0;
    if (x02[m] < x012[m]) {
        var [x_left, x_right] = [x02, x012];
        var [iz_left, iz_right] = [z02, z012];
      } else {
        var [x_left, x_right] = [x012, x02];
        var [iz_left, iz_right] = [z012, z02];
      }
    
      // Draw horizontal segments.
      for (let y = p0.y; y <= p2.y; y++) {
        let [xl, xr] = [x_left[y - p0.y] | 0, x_right[y - p0.y] | 0];
    
        // Interpolate attributes for this scanline.
        let [zl, zr] = [iz_left[y - p0.y], iz_right[y - p0.y]];
        let zscan = interpolate(xl, zl, xr, zr);
    
        for (let x = xl; x <= xr; x++) {
          if (UpdateDepthBufferIfCloser(x, y, zscan[x - xl])) {
            putPixel(x, y, {r: triangle.color.r * intensity, g: triangle.color.g * intensity, b: triangle.color.b * intensity});
          }
        }
      }
    // drawWireframeTriangle(
    //     projected[triangle.v[0]], 
    //     projected[triangle.v[1]], 
    //     projected[triangle.v[2]], 
    //     {r: triangle.color.r * 0.75, g: triangle.color.g * 0.75, b: triangle.color.b * 0.75});
}




const Identity4x4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


function makeOYRotationMatrix(rotation) {
    const degree = rotation.y;
    const cos = Math.cos(degree*Math.PI/180.0);
    const sin = Math.sin(degree*Math.PI/180.0);
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
        [0, 0, 0,             1]
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
    let result = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];

    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            for (let k = 0; k < 4; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

// multiply matrix and 4d vector
function multiplyMatrixVector(matrix, vector) {
    let result = [0, 0, 0, 0];
    let vec = [vector.x, vector.y, vector.z, vector.w];

    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            result[i] += matrix[i][j] * vec[j];
        }
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


const signedDistance = (vertex, plane) => {
    const normal = plane.normal;
    const d = plane.d;
    return normal.x * vertex.x + normal.y * vertex.y + normal.z * vertex.z + d;
}

// send 
const clipTrianglesAgainstPlane = (model, plane, vertices) => {
    
    let clipped_triangles = [];

    for (let i = 0; i < model.triangles.length; i++) {
        let triangle_new = model.triangles[i];
        let clipped_triangle = clipTriangle(triangle_new, vertices, plane);

        if (clipped_triangle === null) {
            continue; // Skip this triangle instead of returning null
        }
        clipped_triangles.push(clipped_triangle);
    }
    return clipped_triangles.length > 0 ? clipped_triangles : null;
}


const clipTriangle = (triangle, vertices, plane) => {
    // console.log(vertices[triangle.v[0]]);

    // Getting signed distance for each vertex of the triangle with the plane

    let d0 = signedDistance(vertices[triangle.v[0]], plane);
    let d1 = signedDistance(vertices[triangle.v[1]], plane) ;
    let d2 = signedDistance(vertices[triangle.v[2]], plane);
 
    let in0 = d0 > 0;
    let in1 = d1 > 0;
    let in2 = d2 > 0;
    
    // Count how many vertices are in front of the plane
    let in_count = (in0 ? 1 : 0) + (in1 ? 1 : 0) + (in2 ? 1 : 0);


    if (in_count == 0) {
        // Nothing to do - the triangle is fully clipped out.
      } else if (in_count == 3) {
        // The triangle is fully in front of the plane.
        return triangle;
    } else if (in_count == 1) {
        // The triangle has one vertex in. Output is one clipped triangle.
      } else if (in_count == 2) {
        // The triangle has two vertices in. Output is two clipped triangles.
      }    
}




const clipInstance = (model, planes, transformMatrix, scale) => {
    // Transform center and radius of the bounding sphere

    let center = multiplyMatrixVector(transformMatrix, new Vec4(model.bounding_center.x, model.bounding_center.y, model.bounding_center.z, 1));
    let radius = model.bounding_radius * scale.x;

    for (let i = 0; i < planes.length; i++) {
        let plane = planes[i];

        let distance = signedDistance(center, plane);

        // let clipped_instance = clipInstanceAgainstPlane(model, plane, transformMatrix);
        if (distance < -radius) {
            // If the instance is completely outside the plane return null
            return null;
        }
    }

    // applying model view transform
    let vertices = []
    for (const vertex of model.vertices) {
        let projectedVertex = multiplyMatrixVector(transformMatrix, new Vec4(vertex.x, vertex.y, vertex.z, 1))
        vertices.push(projectedVertex)
    }

    // copy of triangles to avoid modifying the original
    let triangles = [...model.triangles];
    // Clip triangles against planes
    for (let plane of planes) {
        let clipped_triangles = clipTrianglesAgainstPlane(model, plane, vertices);
        if (clipped_triangles === null) {
            return null;
        }
        triangles = clipped_triangles;
    }

    return new Model(
        model.name,
        vertices,
        triangles,
        center,
        radius
    );

}

const makeCameraMatrix = (position, rotation) => {
    const translationMatrix = makeTranslationMatrix(position);
    const rotationMatrix = transposeMatrix(makeOYRotationMatrix(rotation));
    const cameraMatrix = multiplyMatrices(rotationMatrix, translationMatrix);
    return cameraMatrix;
}

const renderObject = (instance) => {
    let projected = []
    const instancedModel = instance;
    // console.log(instancedModel)
    
    for (const vertex of instancedModel.vertices) {
        
        let projectedVertex = projectVertex(vertex);
        // console.log(vertex, projectedVertex)
        // console.log(multiplyMatrixVector(vertex))
        projected.push(projectedVertex)

    }


    for (const triangle of instancedModel.triangles) {
        // console.log(triangle)
        renderTriangle(triangle, projected, instancedModel.vertices);
    }
}

const renderScene = () => {
    const cameraMatrix = makeCameraMatrix(new Vertex(-camera.position.x, -camera.position.y, -camera.position.z), camera.rotation);
    // console.log(cameraMatrix);
    for (const cube of cubeInstance) {
        let cubeRotationMatrix = cube.transform.rotation.y !== 0 ? makeOYRotationMatrix(cube.transform.rotation) : Identity4x4
        let cubeTransformMatrix = multiplyMatrices(
            makeTranslationMatrix(cube.transform.translation), 
            multiplyMatrices( 
                cubeRotationMatrix,
                makeScaleMatrix(cube.transform.scale)
            )
        );
        let transformMatrix = multiplyMatrices(cameraMatrix, cubeTransformMatrix);
        let clipped_cube_model = clipInstance(cube.model, clipping_planes, transformMatrix, cube.transform.scale);
        if (clipped_cube_model != null) {
            renderObject(clipped_cube_model);
        }
    }
}

const v0 = new Vertex(1, 1, 1)
const v1 = new Vertex(-1, 1, 1)
const v2 = new Vertex(-1, -1, 1)
const v3 = new Vertex(1, -1, 1)

const v4 = new Vertex(1, 1, -1)
const v5 = new Vertex(-1, 1, -1)
const v6 = new Vertex(-1, -1, -1)
const v7 = new Vertex(1, -1, -1)

const vertices = [v0, v1, v2, v3, v4, v5, v6, v7]

// A Light.
const LT_AMBIENT = 0;
const LT_POINT = 1;
const LT_DIRECTIONAL = 2;

const LM_DIFFUSE = 1;
const LM_SPECULAR = 2;

const SM_FLAT = 0;
const SM_GOURAUD = 1;
const SM_PHONG = 2;

const lights = [
    new Light(LT_AMBIENT, 0.2, new Vertex(0, 0, 0)),
    new Light(LT_DIRECTIONAL, 0.2, new Vertex(-1, 0, 1)),
    new Light(LT_POINT, 0.6, new Vertex(-3, 2, -10))
];
  

// Triangle index for cube
const triangle_index = [
    new Triangle([0, 1, 2], RED),
    new Triangle([0, 2, 3], RED),
    new Triangle([1, 5, 6], YELLOW),
    new Triangle([1, 6, 2], YELLOW),
    new Triangle([2, 6, 7], CYAN),
    new Triangle([2, 7, 3], CYAN),
    new Triangle([4, 0, 3], GREEN),
    new Triangle([4, 1, 0], PURPLE),
    new Triangle([4, 3, 7], GREEN),
    new Triangle([4, 5, 1], PURPLE),
    new Triangle([5, 4, 7], BLUE),
    new Triangle([5, 7, 6], BLUE),
  ];
  

const cubeModel = new Model("cube", vertices, triangle_index, new Vertex(0, 0, 0), Math.sqrt(3));

const cubeA = new Instance(cubeModel, new Transform({x: 0.75,y: 0.75,z: 0.75}, {x: 0, y: 0, z: 0}, {x: -1.5, y: 0, z: 7}))
const cubeB = new Instance(cubeModel, new Transform({x: 1,y: 1,z: 1}, {x: 0, y: 195, z: 0}, {x: 1.25, y: 2.5, z: 7.5}))
const sphereA = new Instance(GenerateSphere(15, GREEN), new Transform({x: 1.5, y: 1.5, z: 1.5}, {x: 0, y: 0, z: 0}, {x: 1.75, y: -0.5, z: 7}))
const cubeInstance = [cubeA, cubeB, sphereA]

let camera = new Camera(new Vertex(-3, 1,2), {x: 0, y: -30, z: 0})

let s2 = 1.0/Math.sqrt(2);

const clipping_planes = [
    new Plane(new Vertex(0,0,1), -1), // Near
    new Plane(new Vertex(s2,0,s2), 0), // left
    new Plane(new Vertex(-s2,0,s2), 0), // right
    new Plane(new Vertex(0,-s2,s2), 0), // top
    new Plane(new Vertex(0,s2,s2), 0), // bottom
]

console.log(depth_buffer)

renderScene()

updateCanvas();