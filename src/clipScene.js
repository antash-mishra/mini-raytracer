const clipScene = (instances, planes) => {
    clipped_instance = [];
    for (let i = 0; i < instances.length; i++) {
        let instance = instances[i];
        let clipped_instance = clipInstance(instance, planes);
        if (clipped_instance !== null) {
            // Check if the clipped instance is not empty
            clipped_instances.push(clipped_instance);
        }
    }
    return clipped_instances;

}

const clipInstance = (instance, planes) => {
    for (let i = 0; i < planes.length; i++) {
        let plane = planes[i];
        let clipped_instance = clipInstanceAgainstPlane(instance, plane);
        if (clipped_instance === null) {
            return null;
        }
        instance = clipped_instance;
    }
    return instance;
}

const clipInstanceAgainstPlane = (instance, plane) => {
    // Getting signed Distance to the plane
    let distance = signedDistance(instance.bounding_center, plane);
    let radius = instance.bounding_radius;

    // If the instance is completely inside the plane return instance
    if (distance > radius) {
        return instance;
    }
    else if (distance < -radius) {
        // If the instance is completely outside the plane return null
        return null;
    }
    else {
        // If the instance is partially inside the plane, clip it
        let clipped_instance = clipTrianglesAgainstPlane(instance.triangles, plane);
        return clipped_instance;
    }
}


const clipTrianglesAgainstPlane(triangles, plane) => {
    let clipped_triangles = [];
    for (let i = 0; i < triangles.length; i++) {
        let triangle = triangles[i];
        let clipped_triangle = clipTriangle(triangle, plane);
        if (clipped_triangle !== null) {
            clipped_triangles.push(clipped_triangle);
        }
    }
    return clipped_triangles;
}

const clipTriangle = (triangle, plane) => {
    // Getting signed distance for each vertex of the triangle with the plane
    let d0 = signedDistance(triangle[0], plane) > 0;
    let d1 = signedDistance(triangle[1], plane) > 0;
    let d2 = signedDistance(triangle[2], plane) > 0;

    let in_count = d0 + d1 + d2;

    if (in_count === 3) {
        // If all vertices are inside the plane, return the triangle
        return triangle;
    }
}