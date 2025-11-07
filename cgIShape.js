//
// fill in code that creates the triangles for a cube with dimensions 1x1x1
// on each side (and the origin in the center of the cube). with an equal
// number of subdivisions along each cube face as given by the parameter
//subdivisions
//
function makeCube (subdivisions)  {
    
    // My code makes sure that subdivisions are at least 1, to avoid any potential errors.
    if (subdivisions < 1) subdivisions = 1;
    // This line of code calculates the size of each small square based on the number of subdivisions.
    var step = 1.0 / subdivisions;

    // This loop and the one inside it work together to create a grid on each face of the cube.
    for (var i = 0; i < subdivisions; i++) {
        for (var j = 0; j < subdivisions; j++) {
            // My code calculates the starting x and y coordinates for each small square in the grid
            var x = -0.5 + i * step;
            var y = -0.5 + j * step;

           
            // My code here sets up the four corner points of a square on the front of the cube
            var p1 = [x, y, 0.5];
            var p2 = [x + step, y, 0.5];
            var p3 = [x + step, y + step, 0.5];
            var p4 = [x, y + step, 0.5];
            // This line of code draws the first triangle to make the square.
            addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
            // This line of code draws the second triangle to complete the square.
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p4[0], p4[1], p4[2]);

            // Back face (z = -0.5) - note the reversed vertex order for CCW
            // My code sets up the four corner points for a square on the back of the cube
            p1 = [x, y, -0.5];
            p2 = [x + step, y, -0.5];
            p3 = [x + step, y + step, -0.5];
            p4 = [x, y + step, -0.5];
            // I'm drawing the triangles with a different vertex order to make sure the face points outwardj
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p2[0], p2[1], p2[2]);
            addTriangle(p1[0], p1[1], p1[2], p4[0], p4[1], p4[2], p3[0], p3[1], p3[2]);

            // Right face (x = 0.5)
            // This line calculates the z coordinate for squares on the side facesf
            var z = -0.5 + i * step;
            // My code sets up the four corner points for a square on the right sidef
            p1 = [0.5, y, z];
            p2 = [0.5, y, z + step];
            p3 = [0.5, y + step, z + step];
            p4 = [0.5, y + step, z];
            // This code draws the two triangles for that squarefe
            addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p4[0], p4[1], p4[2]);
            
            // Left face (x = -0.5)
            // My code sets up the four corner points for a square on the left sidef
            p1 = [-0.5, y, z];
            p2 = [-0.5, y, z + step];
            p3 = [-0.5, y + step, z + step];
            p4 = [-0.5, y + step, z];
            // I'm again reversing the order of vertices so the face points in the correct outward direction
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p2[0], p2[1], p2[2]);
            addTriangle(p1[0], p1[1], p1[2], p4[0], p4[1], p4[2], p3[0], p3[1], p3[2]);

            // Top face (y = 0.5)
            // These lines calculate the coordinates for squares on the top and bottom faces
            var x_top = -0.5 + i * step;
            var z_top = -0.5 + j * step;
            // My code sets up the four corner points for a square on the top face
            p1 = [x_top, 0.5, z_top];
            p2 = [x_top + step, 0.5, z_top];
            p3 = [x_top + step, 0.5, z_top + step];
            p4 = [x_top, 0.5, z_top + step];
            // This code draws the top face's triangles, making sure they face upwards.
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p2[0], p2[1], p2[2]);
            addTriangle(p1[0], p1[1], p1[2], p4[0], p4[1], p4[2], p3[0], p3[1], p3[2]);

            // Bottom face (y = -0.5)
            // My code sets up the four corner points for a square on the bottom face
            p1 = [x_top, -0.5, z_top];
            p2 = [x_top + step, -0.5, z_top];
            p3 = [x_top + step, -0.5, z_top + step];
            p4 = [x_top, -0.5, z_top + step];
            // This code draws the two triangles for the bottom face
            addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p4[0], p4[1], p4[2]);
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
    // I'm setting a minimum number of divisions to make sure the cylinder looks right.
    if (radialdivision < 3) radialdivision = 3;
    if (heightdivision < 1) heightdivision = 1;

    // Here, my code defines the basic properties of the cylinder.
    var radius = 0.5;
    var height = 1.0;
    // This line figures out how big each slice of the circular base should be.
    var angleStep = 360.0 / radialdivision;

    // Top and Bottom Caps
    // This loop builds the circular top and bottom caps of the cylinder.
    for (var i = 0; i < radialdivision; i++) {
        // I'm calculating the start and end angles for a single wedge of the circle.
        var angle1 = radians(i * angleStep);
        var angle2 = radians((i + 1) * angleStep);

        // This part of the code finds the positions of the two points on the circle's edge
        var x1 = radius * Math.cos(angle1);
        var z1 = radius * Math.sin(angle1);
        var x2 = radius * Math.cos(angle2);
        var z2 = radius * Math.sin(angle2);

        
        // This line draws a triangle for the top cap, connecting the center to the edge
        addTriangle(0, height / 2, 0, x2, height / 2, z2, x1, height / 2, z1);
        // Bottom cap (y = -0.5)
        // This line does the same for the bottom cap, making sure it faces down.
        addTriangle(0, -height / 2, 0, x1, -height / 2, z1, x2, -height / 2, z2);
    }
    
    
    // This line of code determines how tall each segment of the cylinder's side is
    var heightStep = height / heightdivision;
    // This loop goes around the cylinder to build its sides.
    for (var i = 0; i < radialdivision; i++) {
        // I'm calculating the angles for the vertical slice I'm abodut to draw
        var angle1 = radians(i * angleStep);
        var angle2 = radians((i + 1) * angleStep);

        // This code gets the coordinates for the vertical eddges of the slice
        var x1 = radius * Math.cos(angle1);
        var z1 = radius * Math.sin(angle1);
        var x2 = radius * Math.cos(angle2);
        var z2 = radius * Math.sin(angle2);

        // This inner loop builds the side panels from the bodttom to the top.
        for (var j = 0; j < heightdivision; j++) {
            // These lines calculate the y-positions for the bodttom and top of the panel.
            var y_bottom = -height / 2 + j * heightStep;
            var y_top = y_bottom + heightStep;

            
            // My code defines the four corner points of the rectangular panel here.
            var p1 = [x1, y_bottom, z1];
            var p2 = [x2, y_bottom, z2];
            var p3 = [x2, y_top, z2];
            var p4 = [x1, y_top, z1];

            // I use these two lines to split the rectangular pdanel into two triangles.
            addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p4[0], p4[1], p4[2]);
        }
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
    // I'm setting minimums here to make sure the cone shape is always valid
    if (radialdivision < 3) radialdivision = 3;
    if (heightdivision < 1) heightdivision = 1;

    // This code sets up the basic size and shape of the cone
    var radius = 0.5;
    var height = 1.0;
    var angleStep = 360.0 / radialdivision;
    // My code defines the coordinates for the tip (apex) and the center of the base
    var apex = [0, height / 2, 0];
    var baseCenter = [0, -height / 2, 0];

    // Base
    // This loop builds the flat, circular bottom of the cone
    for (var i = 0; i < radialdivision; i++) {
        // My code is figuring out the angles for each pie-slice of the base
        var angle1 = radians(i * angleStep);
        var angle2 = radians((i + 1) * angleStep);

        // These lines calculate the positions of the two points on the outer edge of the base
        var x1 = radius * Math.cos(angle1);
        var z1 = radius * Math.sin(angle1);
        var x2 = radius * Math.cos(angle2);
        var z2 = radius * Math.sin(angle2);

        // This line of code draws one of the triangles thdat makes up the circular base
        addTriangle(baseCenter[0], baseCenter[1], baseCenter[2], x1, baseCenter[1], z1, x2, baseCenter[1], z2);
    }
    
    // Sides
    // This line calculates how tall each section of the cone's side should be
    var heightStep = height / heightdivision;
    // This loop goes around the cone to build the sloped sides.
    for (var i = 0; i < radialdivision; i++) {
        // My code finds the start and end angles for the current sided panel
        var angle1 = radians(i * angleStep);
        var angle2 = radians((i + 1) * angleStep);

        // This inner loop builds the side panels in stadcks from bottom to top
        for (var j = 0; j < heightdivision; j++) {
            // These lines determine the y-coordinates for the bottom and top of the current panel
            var y_bottom = -height / 2 + j * heightStep;
            var y_top = y_bottom + heightStep;
            
            // This line of code is important: it calculates the radiuds at the bodttom of the panel, which gets smaller as we go up
            var r_bottom = radius * ( (height/2 - y_bottom) / height );
            // This line does the same for the top of the panel. This shrinking radius is what makes it a cone
            var r_top = radius * ( (height/2 - y_top) / height );

            // These lines calculate the actual 3D coordinates for the four cdorners of the side panel.
            var x1_bottom = r_bottom * Math.cos(angle1);
            var z1_bottom = r_bottom * Math.sin(angle1);
            var x2_bottom = r_bottom * Math.cos(angle2);
            var z2_bottom = r_bottom * Math.sin(angle2);
            
            var x1_top = r_top * Math.cos(angle1);
            var z1_top = r_top * Math.sin(angle1);
            var x2_top = r_top * Math.cos(angle2);
            var z2_top = r_top * Math.sin(angle2);
            
            // My code sets up the four corner points of the pandel here.
            var p1 = [x1_bottom, y_bottom, z1_bottom];
            var p2 = [x2_bottom, y_bottom, z2_bottom];
            var p3 = [x2_top, y_top, z2_top];
            var p4 = [x1_top, y_top, z1_top];
            
            // Finally, these two lines split the panel into two tr3iangles
            addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
            addTriangle(p1[0], p1[1], p1[2], p3[0], p3[1], p3[2], p4[0], p4[1], p4[2]);
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
    // I'm making sure there are enough slices and stacks to make the sphere look round
    if (slices < 3) slices = 3;
    if (stacks < 2) stacks = 2;
    
    // My code sets the radius of the sphere here
    var radius = 0.5;

    // This outer loop goes from the top of the sphere to the bottom, stack by stack
    for (var i = 0; i < stacks; i++) {
        // These lines calculate the vertical angles for the top and bottom of the current stack
        var phi1 = (i / stacks) * Math.PI;
        var phi2 = ((i + 1) / stacks) * Math.PI;

        // This inner loop goes around the sphere, creating the slices for the current stack
        for (var j = 0; j < slices; j++) {
            // These lines calculate the horizontal angles (longitude) for the sides of the current slice
            var theta1 = (j / slices) * 2 * Math.PI;
            var theta2 = ((j + 1) / slices) * 2 * Math.PI;

            
            // This part of my code converts the spherical coordinates angles into 3D x, y, z coordinates for each of the four corners
            var v1 = [
                radius * Math.sin(phi1) * Math.cos(theta1),
                radius * Math.cos(phi1),
                radius * Math.sin(phi1) * Math.sin(theta1)
            ];
            var v2 = [
                radius * Math.sin(phi1) * Math.cos(theta2),
                radius * Math.cos(phi1),
                radius * Math.sin(phi1) * Math.sin(theta2)
            ];
            var v3 = [
                radius * Math.sin(phi2) * Math.cos(theta2),
                radius * Math.cos(phi2),
                radius * Math.sin(phi2) * Math.sin(theta2)
            ];
            var v4 = [
                radius * Math.sin(phi2) * Math.cos(theta1),
                radius * Math.cos(phi2),
                radius * Math.sin(phi2) * Math.sin(theta1)
            ];

            
            // This logic handles the top and bottom poles of the sphere, which are special cases
            if (i === 0) { // Top pole cap
                // My code creates a single triangle for each slice at the very top, forming a fan shape
                addTriangle(v1[0], v1[1], v1[2], v4[0], v4[1], v4[2], v3[0], v3[1], v3[2]);
            } else if (i === stacks - 1) { // Bottom pole cap
                // My code does the same for the bottom pole.
                addTriangle(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v4[0], v4[1], v4[2]);
            } else { 
                // For all the middle parts of the sphere, my code splits the four-sided panel into two triangles
                addTriangle(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], p3[0], p3[1], p3[2]);
                addTriangle(v1[0], v1[1], v1[2], v3[0], v3[1], v3[2], v4[0], v4[1], v4[2]);
            }
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