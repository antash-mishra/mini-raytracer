class Sphere {
    constructor(center, radius, color, specular) {
        this.center = center;
        this.radius = radius;
        this.color = color;
        this.specular = specular;
    }
}

class Light {
    constructor(type, intensity, position, direction) {
        this.type = type;
        this.intensity = intensity;
        this.position = position;
        this.direction = direction;
    }
}

const sizes = {
    width: window.innerWidth,
    height: window.innerHeight,
    originX: 0,
    originY: 0,
    originZ: 0,
    viewportWidth: 1,
    viewportHeight: 1,
    projectionPlaneD: 1,
}

const O = {
    x: sizes.originX,
    y: sizes.originY,
    z: sizes.originZ,
}

const BACKGROUND_COLOR = { r: 255, g: 255, b: 255 }

const sphereA = new Sphere({ x: 0, y: -1, z: 3 }, 1, { r: 255, g: 0, b: 0 }, 500)
const sphereB = new Sphere({ x: 2, y: 0, z: 4 }, 1, { r: 0, g: 0, b: 255 }, 500)
const sphereC = new Sphere({ x: -2, y: 0, z: 4 }, 1, { r: 0, g: 255, b: 0 }, 10)
const sphereD = new Sphere({ x: 0, y: -5001, z: 0 }, 5000, { r: 255, g: 255, b: 0 }, 1000)

const spheres = [sphereA, sphereB, sphereC, sphereD]

// Light 
const ambient  = new Light("ambient", 0.2, null, null)
const point = new Light("point", 0.6, {x: 2, y: 1, z: 0}, null)
const directional = new Light("directional", 0.2, null, {x: 1, y: 4, z: 4})

const lights = [ambient, point, directional]

const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
let canvas_buffer = ctx.getImageData(0, 0, sizes.width, sizes.height);
ctx.fillStyle = "green";

canvas.width = sizes.width;
canvas.height = sizes.height;

const canvasToViewport = (a, b) => {
    return { x: a * (sizes.viewportWidth / sizes.width), y: b * (sizes.viewportHeight / sizes.height), z: sizes.projectionPlaneD };
}

const dot = (a, b) => {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

const magnitude = (v) => {
    const result = Math.sqrt(dot(v, v))
    return result;
}

function intersectRaySphere(O, D, sphere) {
    const radius = sphere.radius;
    const centerToOrigin = {
        x: O.x - sphere.center.x,
        y: O.y - sphere.center.y,
        z: O.z - sphere.center.z
    };

    const a = dot(D, D);
    const b = 2 * dot(D, centerToOrigin);
    const c = dot(centerToOrigin, centerToOrigin) - radius * radius;
    const discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return [Infinity, Infinity];
    }

    const t1 = (-b + Math.sqrt(discriminant)) / (2 * a);
    const t2 = (-b - Math.sqrt(discriminant)) / (2 * a);
    return [t1, t2];
}

// Lightning Calculation
// P => Position
// N => Normal
// V => Point of View by user
// specular => specular Intensity
const computeLighting = (P, N, V, specular) => {
    let i = 0.0;
    for (const light of lights) {
        if (light.type === "ambient") {
            i += light.intensity;
        } else {
            let L = { x: 0, y: 0, z: 0 };
            if (light.type === "point") {
                L = {
                    x: light.position.x - P.x,
                    y: light.position.y - P.y,
                    z: light.position.z - P.z
                };
            } else {
                L = light.direction
            }
            // Diffuse Light Equation
            const n_dot_l = dot(N, L);
            if (n_dot_l > 0) {
                i += light.intensity * n_dot_l / (magnitude(N) * magnitude(L));
            }

            // Specular Light Eq
            if (specular != -1) {
                const reflection = { x: 2 * n_dot_l * N.x - L.x, y: 2 * n_dot_l * N.y - L.y, z: 2 * n_dot_l * N.z - L.z };
                const r_dot_v = dot(reflection, V);
                if (r_dot_v > 0) {
                    i += light.intensity * Math.pow(r_dot_v / (magnitude(reflection) * magnitude(V)), specular);
                }
            }
            
        }
    }
    return i;
}

const traceRay = (O, D, t_min, t_max) => {
    let closest_t = Infinity;
    let closest_sphere = null;
    for (const sphere of spheres) {

        const [t1, t2] = intersectRaySphere(O, D, sphere);
        
        if (t1 > t_min && t1 < t_max && t1 < closest_t) {
            closest_t = t1;
            closest_sphere = sphere;
        }

        if (t2 > t_min && t2 < t_max && t2 < closest_t) {
            closest_t = t2;
            closest_sphere = sphere;
        }
    }
    if (closest_sphere == null) {
        return BACKGROUND_COLOR;
    }

    // return closest_sphere.color;
    // Position
    if (closest_t === Infinity) {
        return BACKGROUND_COLOR;
    }

    let P = {x:  O.x + closest_t * D.x, y: O.y + closest_t * D.y, z: O.z + closest_t * D.z};
    
    // Normal
    let N = {x: P.x - closest_sphere.center.x, y: P.y - closest_sphere.center.y, z: P.z - closest_sphere.center.z}
    N = {x: N.x / magnitude(N), y: N.y / magnitude(N), z: N.z / magnitude(N)}

    // Lighting
    const lighting = computeLighting(P, N, {x: -D.x, y: -D.y, z: -D.z}, closest_sphere.specular)
    return { r: closest_sphere.color.r * lighting, g: closest_sphere.color.g * lighting, b: closest_sphere.color.b * lighting }
}

// const putPixel = (x, y, color) => {
//     x += sizes.width / 2;
//     y = sizes.height / 2 - y;
//     const imageData = ctx.createImageData(1, 1);
    
//     const data = imageData.data;
    
    
//     data[0] = color.r;
//     data[1] = color.g;
//     data[2] = color.b;
//     data[3] = 255;
//     ctx.putImageData(imageData, x, y);
// }



const putPixel = (x, y, color) => {
    x += sizes.width / 2;
    y = sizes.height / 2 - y-1;
    if (x < 0 || x >= sizes.width || y < 0 || y >= sizes.height) {
        return;
    }
    const imageData = canvas_buffer
    const data = imageData.data;
    let offset = 4*(x + canvas_buffer.width*y);
    data[offset++] = color.r;
    data[offset++] = color.g;
    data[offset++] = color.b;
    data[offset++] = 255;
}

function updateCanvas() {
    ctx.putImageData(canvas_buffer, 0, 0);
}

// Fix the rendering loop to use proper coordinates
for (let x = -sizes.width/2; x <= sizes.width/2; x++) {
    for (let y = -sizes.height/2; y <= sizes.height/2; y++) {
        // Convert from screen space to viewport space
        let D = canvasToViewport(x, y);
        let color = traceRay(O, D, 1, Infinity);
        // Ensure proper rendering position
        putPixel(x, y, color);
    }
}

updateCanvas();