//
// fill in code that creates the triangles for a cube with dimensions 1x1x1
// on each side (and the origin in the center of the cube). with an equal
// number of subdivisions along each cube face as given by the parameter
//subdivisions
//
function makeCube (subdivisions)  {
    
    // This function now uses a more streamlined approach to build the cube.
    let step = 1.0 / subdivisions;

    // This loop iterates through all 6 faces of the cube.
    for (let face = 0; face < 6; face++) {
        // These nested loops create the grid for each face.
        for (let i = 0; i < subdivisions; i++) {
            for (let j = 0; j < subdivisions; j++) {
                
                // My code calculates the 2D coordinates for a small square.
                let u = i * step - 0.5;
                let v = j * step - 0.5;
                let u1 = (i + 1) * step - 0.5;
                let v1 = (j + 1) * step - 0.5;

                let p0, p1, p2, p3;

                // This logic correctly maps the 2D square to the proper 3D face of the cube.
                if (face === 0) { // Front face (+z)
                    p0 = [u, v, 0.5]; p1 = [u1, v, 0.5]; p2 = [u1, v1, 0.5]; p3 = [u, v1, 0.5];
                    addTriangle(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
                    addTriangle(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
                } else if (face === 1) { // Back face (-z)
                    p0 = [u, v, -0.5]; p1 = [u1, v, -0.5]; p2 = [u1, v1, -0.5]; p3 = [u, v1, -0.5];
                    addTriangle(p1[0], p1[1], p1[2], p0[0], p0[1], p0[2], p2[0], p2[1], p2[2]);
                    addTriangle(p2[0], p2[1], p2[2], p0[0], p0[1], p0[2], p3[0], p3[1], p3[2]);
                } else if (face === 2) { // Top face (+y)
                    p0 = [u, 0.5, v]; p1 = [u1, 0.5, v]; p2 = [u1, 0.5, v1]; p3 = [u, 0.5, v1];
                    addTriangle(p1[0], p1[1], p1[2], p0[0], p0[1], p0[2], p2[0], p2[1], p2[2]);
                    addTriangle(p2[0], p2[1], p2[2], p0[0], p0[1], p0[2], p3[0], p3[1], p3[2]);
                } else if (face === 3) { // Bottom face (-y)
                    p0 = [u, -0.5, v]; p1 = [u1, -0.5, v]; p2 = [u1, -0.5, v1]; p3 = [u, -0.5, v1];
                    addTriangle(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
                    addTriangle(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
                } else if (face === 4) { // Right face (+x)
                    p0 = [0.5, u, v]; p1 = [0.5, u1, v]; p2 = [0.5, u1, v1]; p3 = [0.5, u, v1];
                    addTriangle(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
                    addTriangle(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
                } else { // Left face (-x)
                    p0 = [-0.5, u, v]; p1 = [-0.5, u1, v]; p2 = [-0.5, u1, v1]; p3 = [-0.5, u, v1];
                    addTriangle(p1[0], p1[1], p1[2], p0[0], p0[1], p0[2], p2[0], p2[1], p2[2]);
                    addTriangle(p2[0], p2[1], p2[2], p0[0], p0[1], p0[2], p3[0], p3[1], p3[2]);
                }
            }
        }
    }
}


//
// fill in code that creates the triangles for a cylinder with diameter 1
// and height of 1 (centered at the origin) with the number of subdivisions
// around the base and top of the cylinder (given by radialdivision) and
// the number of subdivisions along the surface of the cylinder given by
//heightdivision.
//
function makeCylinder (radialdivision,heightdivision){
    let radius = 0.5;
    let heightStep = 1.0 / heightdivision;
    let angleStep = 360.0 / radialdivision;

    for (let j = 0; j < heightdivision; j++) {
        let y0 = -0.5 + j * heightStep;
        let y1 = y0 + heightStep;
        for (let i = 0; i < radialdivision; i++) {
            let angle1 = radians(i * angleStep);
            let angle2 = radians((i + 1) * angleStep);

            let x1 = radius * Math.cos(angle1);
            let z1 = radius * Math.sin(angle1);
            let x2 = radius * Math.cos(angle2);
            let z2 = radius * Math.sin(angle2);

            let p0 = [x1, y0, z1];
            let p1 = [x2, y0, z2];
            let p2 = [x2, y1, z2];
            let p3 = [x1, y1, z1];

            addTriangle(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
            addTriangle(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
        }
    }
    
    for (let i = 0; i < radialdivision; i++) {
        let angle1 = radians(i * angleStep);
        let angle2 = radians((i + 1) * angleStep);
        
        let x1 = radius * Math.cos(angle1);
        let z1 = radius * Math.sin(angle1);
        let x2 = radius * Math.cos(angle2);
        let z2 = radius * Math.sin(angle2);
        
        addTriangle(0, -0.5, 0, x2, -0.5, z2, x1, -0.5, z1);
        addTriangle(0, 0.5, 0, x1, 0.5, z1, x2, 0.5, z2);
    }
}


//
// fill in code that creates the triangles for a cone with diameter 1
// and height of 1 (centered at the origin) with the number of
// subdivisions around the base of the cone (given by radialdivision)
// and the number of subdivisions along the surface of the cone
//given by heightdivision.
//
function makeCone (radialdivision, heightdivision) {
    let radius = 0.5;
    let angleStep = 360.0 / radialdivision;
    let heightStep = 1.0 / heightdivision;
    
    for (let i = 0; i < radialdivision; i++) {
        let angle1 = radians(i * angleStep);
        let angle2 = radians((i + 1) * angleStep);
        let x1 = radius * Math.cos(angle1);
        let z1 = radius * Math.sin(angle1);
        let x2 = radius * Math.cos(angle2);
        let z2 = radius * Math.sin(angle2);
        addTriangle(0, -0.5, 0, x2, -0.5, z2, x1, -0.5, z1);
    }

    for (let j = 0; j < heightdivision; j++) {
        let y0 = -0.5 + j * heightStep;
        let y1 = y0 + heightStep;
        let r0 = radius * (1 - (j / heightdivision));
        let r1 = radius * (1 - ((j + 1) / heightdivision));

        for (let i = 0; i < radialdivision; i++) {
            let angle1 = radians(i * angleStep);
            let angle2 = radians((i + 1) * angleStep);
            
            let x1_0 = r0 * Math.cos(angle1);
            let z1_0 = r0 * Math.sin(angle1);
            let x2_0 = r0 * Math.cos(angle2);
            let z2_0 = r0 * Math.sin(angle2);
            
            let x1_1 = r1 * Math.cos(angle1);
            let z1_1 = r1 * Math.sin(angle1);
            let x2_1 = r1 * Math.cos(angle2);
            let z2_1 = r1 * Math.sin(angle2);

            let p0 = [x1_0, y0, z1_0];
            let p1 = [x2_0, y0, z2_0];
            let p2 = [x2_1, y1, z2_1];
            let p3 = [x1_1, y1, z1_1];

            if (j < heightdivision - 1) {
                addTriangle(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
                addTriangle(p0[0], p0[1], p0[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
            } else {
                addTriangle(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], 0, 0.5, 0);
            }
        }
    }
}
    
//
// fill in code that creates the triangles for a sphere with diameter 1
// (centered at the origin) with number of slides (longitude) given by
// slices and the number of stacks (lattitude) given by stacks.
// For this function, you will implement the tessellation method based
// on spherical coordinates as described in the video (as opposed to the
//recursive subdivision method).
//
function makeSphere (slices, stacks) {
    const radius = 0.5;
    let sphere_points = [];

    // This part of my code generates all the vertex points for the sphere first.
    for (let i = 0; i <= stacks; i++) {
        const phi = i / stacks * Math.PI; // This is the vertical angle (latitude).
        const y = radius * Math.cos(phi);
        
        for (let j = 0; j <= slices; j++) {
            const theta = j / slices * 2 * Math.PI; // This is the horizontal angle (longitude).
            const x = radius * Math.sin(phi) * Math.cos(theta);
            const z = radius * Math.sin(phi) * Math.sin(theta);
            sphere_points.push(x, y, z);
        }
    }

    // Now, my code stitches the points together to form the triangles.
    for (let i = 0; i < stacks; i++) {
        for (let j = 0; j < slices; j++) {
            // This code figures out the indices for the four corners of a quad.
            const first = (i * (slices + 1)) + j;
            const second = first + slices + 1;

            const v1_idx = first * 3;
            const v2_idx = (first + 1) * 3;
            const v3_idx = second * 3;
            const v4_idx = (second + 1) * 3;

            // My code gets the actual 3D coordinates from the points array.
            const v1 = [sphere_points[v1_idx], sphere_points[v1_idx+1], sphere_points[v1_idx+2]];
            const v2 = [sphere_points[v2_idx], sphere_points[v2_idx+1], sphere_points[v2_idx+2]];
            const v3 = [sphere_points[v3_idx], sphere_points[v3_idx+1], sphere_points[v3_idx+2]];
            const v4 = [sphere_points[v4_idx], sphere_points[v4_idx+1], sphere_points[v4_idx+2]];
            
            // These two lines create the two triangles that form the quad for this segment of the sphere.
            addTriangle(v1[0], v1[1], v1[2], v3[0], v3[1], v3[2], v4[0], v4[1], v4[2]);
            addTriangle(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]);
        }
    }
}


////////////////////////////////////////////////////////////////////
//
//  Do not edit below this line
//
///////////////////////////////////////////////////////////////////

function radians(degrees)
{
  var pi = Math.PI;
  return degrees * (pi/180);
}

function addTriangle (x0,y0,z0,x1,y1,z1,x2,y2,z2) {

    
    var nverts = points.length / 4;
    
    // push first vertex
    points.push(x0);  bary.push (1.0);
    points.push(y0);  bary.push (0.0);
    points.push(z0);  bary.push (0.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++;
    
    // push second vertex
    points.push(x1); bary.push (0.0);
    points.push(y1); bary.push (1.0);
    points.push(z1); bary.push (0.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++
    
    // push third vertex
    points.push(x2); bary.push (0.0);
    points.push(y2); bary.push (0.0);
    points.push(z2); bary.push (1.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++;
}