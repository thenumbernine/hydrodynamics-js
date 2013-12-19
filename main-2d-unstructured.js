/*
1D ADM
2D: Burger at least
2D ADM
2D unstructured

sources:
http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/
http://www.cfdbooks.com/cfdcodes.html
"Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro
http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf also does not
*/

var panel;
var canvas;

var xmin = -.5;
var xmax = .5; 
var ymin = -.5;
var ymax = .5;

var useNoise = true;
var noiseAmplitude = .01;

var useCFL = true;
var fixedDT = .2;
var gaussSeidelIterations = 20;

var mouse;

var externalForceX = 0;
var externalForceY = 0;

//interface directions
var dirs = [[1,0], [0,1]];

/*
output:
matrix
eigenvalues
eigenvectors

input:
velocity
hTotal
speedOfSound
*/
function buildEigenstate(offset, matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocityX, velocityY, hTotal, gamma, normalX, normalY) {

	if ((hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)) < 0) {
		console.log('sqrt error');
	}

	//calculate matrix & eigenvalues & vectors at interface from state at interface
	var speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)));
	var tangentX = -normalY;
	var tangentY = normalX;
	var velocityN = velocityX * normalX + velocityY * normalY;
	var velocityT = velocityX * tangentX + velocityY * tangentY;
	var velocitySq = velocityX * velocityX + velocityY * velocityY;	
	
	//eigenvalues: min, mid, max
	eigenvalues[0 + 4 * offset] = velocityN - speedOfSound;
	eigenvalues[1 + 4 * offset] = velocityN;
	eigenvalues[2 + 4 * offset] = velocityN;
	eigenvalues[3 + 4 * offset] = velocityN + speedOfSound;

	//I'm going with http://people.nas.nasa.gov/~pulliam/Classes/New_notes/euler_notes.pdf

	//min eigenvector
	eigenvectors[0 + 4 * (0 + 4 * offset)] = 1;
	eigenvectors[1 + 4 * (0 + 4 * offset)] = velocityX - speedOfSound * normalX;
	eigenvectors[2 + 4 * (0 + 4 * offset)] = velocityY - speedOfSound * normalY;
	eigenvectors[3 + 4 * (0 + 4 * offset)] = hTotal - speedOfSound * velocityN;
	//mid eigenvector (normal)
	eigenvectors[0 + 4 * (1 + 4 * offset)] = 1;
	eigenvectors[1 + 4 * (1 + 4 * offset)] = velocityX;
	eigenvectors[2 + 4 * (1 + 4 * offset)] = velocityY;
	eigenvectors[3 + 4 * (1 + 4 * offset)] = .5 * velocitySq;
	//mid eigenvector (tangent)
	eigenvectors[0 + 4 * (2 + 4 * offset)] = 0;
	eigenvectors[1 + 4 * (2 + 4 * offset)] = tangentX;
	eigenvectors[2 + 4 * (2 + 4 * offset)] = tangentY;
	eigenvectors[3 + 4 * (2 + 4 * offset)] = velocityT;
	//max eigenvector
	eigenvectors[0 + 4 * (3 + 4 * offset)] = 1;
	eigenvectors[1 + 4 * (3 + 4 * offset)] = velocityX + speedOfSound * normalX;
	eigenvectors[2 + 4 * (3 + 4 * offset)] = velocityY + speedOfSound * normalY;
	eigenvectors[3 + 4 * (3 + 4 * offset)] = hTotal + speedOfSound * velocityN;
	
	//calculate eigenvector inverses ... 
	//min row
	eigenvectorsInverse[0 + 4 * (0 + 4 * offset)] = (.5 * (gamma - 1) * velocitySq + speedOfSound * velocityN) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[0 + 4 * (1 + 4 * offset)] = -(normalX * speedOfSound + (gamma - 1) * velocityX) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[0 + 4 * (2 + 4 * offset)] = -(normalY * speedOfSound + (gamma - 1) * velocityY) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[0 + 4 * (3 + 4 * offset)] = (gamma - 1) / (2 * speedOfSound * speedOfSound);
	//mid normal row
	eigenvectorsInverse[1 + 4 * (0 + 4 * offset)] = 1 - .5 * (gamma - 1) * velocitySq / (speedOfSound * speedOfSound);
	eigenvectorsInverse[1 + 4 * (1 + 4 * offset)] = (gamma - 1) * velocityX / (speedOfSound * speedOfSound);
	eigenvectorsInverse[1 + 4 * (2 + 4 * offset)] = (gamma - 1) * velocityY / (speedOfSound * speedOfSound);
	eigenvectorsInverse[1 + 4 * (3 + 4 * offset)] = -(gamma - 1) / (speedOfSound * speedOfSound);
	//mid tangent row
	eigenvectorsInverse[2 + 4 * (0 + 4 * offset)] = -velocityT; 
	eigenvectorsInverse[2 + 4 * (1 + 4 * offset)] = tangentX;
	eigenvectorsInverse[2 + 4 * (2 + 4 * offset)] = tangentY;
	eigenvectorsInverse[2 + 4 * (3 + 4 * offset)] = 0;
	//max row
	eigenvectorsInverse[3 + 4 * (0 + 4 * offset)] = (.5 * (gamma - 1) * velocitySq - speedOfSound * velocityN) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[3 + 4 * (1 + 4 * offset)] = (normalX * speedOfSound - (gamma - 1) * velocityX) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[3 + 4 * (2 + 4 * offset)] = (normalY * speedOfSound - (gamma - 1) * velocityY) / (2 * speedOfSound * speedOfSound);
	eigenvectorsInverse[3 + 4 * (3 + 4 * offset)] = (gamma - 1) / (2 * speedOfSound * speedOfSound);

	//calculate matrix
	var identCheck = [];
	var identBad = false;
	for (var i = 0; i < 4; ++i) {
		for (var j = 0; j < 4; ++j) {
			var s = 0;
			var d = 0;
			for (var k = 0; k < 4; ++k) {
				/** /
				s += eigenvectorsInverse[i + 4 * (k + 4 * offset)] * eigenvectors[k + 4 * (j + 4 * offset)] * eigenvalues[k + 4 * offset];
				identCheck += eigenvectorsInverse[i + 4 * (k + 4 * offset)] * eigenvectors[k + 4 * (j + 4 * offset)];
				/**/
				s += eigenvectors[i + 4 * (k + 4 * offset)] * eigenvalues[k + 4 * offset] * eigenvectorsInverse[k + 4 * (j + 4 * offset)];
				d += eigenvectors[i + 4 * (k + 4 * offset)] * eigenvectorsInverse[k + 4 * (j + 4 * offset)];
				/**/
			}
			matrix[i + 4 * (j + 4 * offset)] = s;
			identCheck[i + 4 * j] = d;
			var epsilon = 1e-5;
			if (Math.abs(d - (i == j ? 1 : 0)) > epsilon) identBad = true;
		}
	}
	if (identBad) {
		console.log('bad eigen basis', identCheck);
	}
	/** /	
	function f32subset(a, o, s) {
		var d = new Float64Array(s);
		for (var i = 0; i < s; ++i) {
			d[i] = a[i+o];
		}
		return d;
	}
	console.log('offset',offset);
	console.log('velocity',velocityX,velocityY);
	console.log('hTotal',hTotal);
	console.log('gamma',gamma);
	console.log('normal',normalX,normalY);

	console.log('eigenvalues:',f32subset(eigenvalues, 4*offset, 4));
	console.log('eigenvectors:',f32subset(eigenvectors, 16*offset, 16));
	console.log('eigenvectors^-1:',f32subset(eigenvectorsInverse, 16*offset, 16));
	console.log('matrix:',f32subset(matrix, 16*offset, 16));
	console.log('e^-1 * e:',identCheck);
	throw 'here';
	/**/
}

var drawToScreenMethods = {
	Density : function() {
		var dataMin = this.state.cells[0].q[0];
		var dataMax = dataMin + 1e-9;
		var lastDataRange = this.lastDataMax - this.lastDataMin;
		for (var i = 0; i < this.state.cells.length; ++i) {	
			var s = this.state.cells[i].q[0];
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Velocity : function() {
		var rho0 = this.state.cells[0].q[0];
		var u0 = this.state.cells[0].q[1] / rho0;
		var v0 = this.state.cells[0].q[2] / rho0;
		var dataMin = Math.sqrt(u0 * u0 + v0 * v0);
		var dataMax = dataMin + 1e-9; 
		var lastDataRange = this.lastDataMax - this.lastDataMin;
		for (var i = 0; i < this.state.cells.length; ++i) {
			var cell = this.state.cells[i];
			var rho = cell.q[0];
			var mx = cell.q[1];
			var my = cell.q[2];
			var s = Math.sqrt(mx * mx + my * my) / rho;
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Energy : function() {
		var q0 = this.state.cells[0].q;
		var dataMin = q0[3] / q0[0]; 
		var dataMax = dataMin + 1e-9; 
		var lastDataRange = this.lastDataMax - this.lastDataMin;
		for (var i = 0; i < this.state.cells.length; ++i) {
			var cell = this.state.cells[i];
			var s = cell.q[3] / cell.q[0];
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Pressure : function() {
		var pressure = this.state.cells[0].pressure;
		var dataMin = pressure[0];
		var dataMax = dataMin + 1e-9;
		var lastDataRange = this.lastDataMax - this.lastDataMin;
		for (var i = 0; i < this.state.cells.length; ++i) {
			var s = this.state.cells[i].pressure;
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	}/*,
	Curl : function() {
		var nx = this.state.nx;
		var q = this.state.q;
		var dataMin = undefined; 
		var dataMax = undefined; 
		var lastDataRange = this.lastDataMax - this.lastDataMin;
		var e = 0;
		var dx = (xmax - xmin) / nx;
		var dy = (ymax - ymin) / nx;
		for (var j = 0; j < nx; ++j) {
			for (var i = 0; i < nx; ++i) {
				var rho = q[0 + 4 * e];
				var mx = q[1 + 4 * e];
				var my = q[2 + 4 * e];
				var mxPrev = this.state.q[1 + 4 * ( ((i+nx-1)%nx) + nx * j )];
				var myPrev = this.state.q[1 + 4 * ( i + nx * ((j+nx-1)%nx) )];
				var dmx = mx - mxPrev;
				var dmy = my - myPrev;
				var s = (dmx / dy - dmy / dx) / rho;
				if (dataMin === undefined || s < dataMin) dataMin = s;
				if (dataMax === undefined || s > dataMax) dataMax = s;
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		return [dataMin, dataMax];
	}*/
};

var integrationMethods = {
	'Forward Euler' : {
		advect : {
			Burgers : function(dt) {
				/**/
				for (var i = 0; i < this.edges.length; ++i) {
					//TODO weight averages by distances from edge ... or volume ... or something
					var edge = this.edges[i];
					if (edge.cells.length == 1) {
						//zero boundaries
						var cell = this.cells[edge.cells[0]];
						edge.u = 0;
					} else {
						var cellA = this.cells[edge.cells[0]];
						var cellB = this.cells[edge.cells[1]];
						var rhoA = cellA.q[0];
						var uA = cellA.q[1] / rhoA;
						var vA = cellA.q[2] / rhoA;
						var rhoB = cellB.q[0];
						var uB = cellB.q[1] / rhoB;
						var vB = cellB.q[2] / rhoB;
						edge.u = .5 * ((uA + uB) * edge.normal[0] + (vA + vB) * edge.normal[1]);
					}
					//if (edge.u !== edge.u) throw 'nan';
				}

				//compute flux and advect for each state vector
				for (var i = 0; i < this.edges.length; ++i) {
					var edge = this.edges[i];
					//TODO slope calc
					for (var j = 0; j < 4; ++j) {
						edge.r[j] = 0;
					}
				}

				//construct flux:
				for (var i = 0; i < this.edges.length; ++i) {
					var edge = this.edges[i];
					if (edge.cells.length == 1) {
						//TODO boundary reflection?
						edge.flux[0] = edge.flux[1] = edge.flux[2] = edge.flux[3] = 0;
					} else {
						var cellA = this.cells[edge.cells[0]];
						var cellB = this.cells[edge.cells[1]];
						if (edge.u >= 0) {	//pull from cell B
							for (var j = 0; j < 4; ++j) {
								edge.flux[j] = edge.u * cellB.q[j];
							}
						} else {
							for (var j = 0; j < 4; ++j) {
								edge.flux[j] = edge.u * cellA.q[j];
							}
						}
						//TODO add flux limiter stuff
					}
				}

				//if (dt !== dt) throw 'nan';
				for (var i = 0; i < this.cells.length; ++i) {
					var cell = this.cells[i];
					for (var j = 0; j < cell.edges.length; ++j) {
						var edge = this.edges[cell.edges[j]];
						var dx = edge.x[0] - cell.x[0];
						var dy = edge.x[1] - cell.x[1];
						var ds = Math.sqrt(dx * dx + dy * dy);
						//if (ds !== ds) throw 'nan';
						if (ds === 0) throw 'divide by zero for cell ' + i + ' and edge ' + cell.edges[j];
						for (var k = 0; k < 4; ++k) {
							//if (cell.q[k] !== cell.q[k]) throw 'nan';
							//if (edge.flux[k] !== edge.flux[k]) throw 'nan';
							cell.q[k] -= dt / ds * edge.flux[k];
							//if (cell.q[k] !== cell.q[k]) throw 'nan';
						}
					}
				}
				/**/
				
				//compute pressure
				for (var i = 0; i < this.cells.length; ++i) {
					var cell = this.cells[i];
					//if (cell.q[0] !== cell.q[0]) throw 'nan';
					//if (cell.q[1] !== cell.q[1]) throw 'nan';
					//if (cell.q[2] !== cell.q[2]) throw 'nan';
					//if (cell.q[0] === 0) throw 'divide by zero error';
					var x = cell.x[0];
					var y = cell.x[1];
					var rho = cell.q[0];
					var u = cell.q[1] / rho;
					var v = cell.q[2] / rho;
					var energyTotal = cell.q[3] / rho;
					var energyKinetic = .5 * (u * u + v * v);
					var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
					var energyThermal = energyTotal - energyKinetic - energyPotential;
					cell.pressure = (this.gamma - 1) * rho * energyThermal;
					//if (cell.pressure !== cell.pressure) throw 'nan';
				}

				//apply external force
				for (var i = 0; i < this.cells.length; ++i) {
					var cell = this.cells[i];
					var rho = cell.q[0];
					cell.q[3] -= dt * (externalForceX * cell.q[1] + externalForceY * cell.q[2]);
					cell.q[1] -= dt * rho * externalForceX;
					cell.q[2] -= dt * rho * externalForceY;
				}

				//apply momentum diffusion = pressure
				/* per-edge method */
				for (var i = 0; i < this.edges.length; ++i) {
					var edge = this.edges[i];
					if (edge.cells.length == 2) {	//one neighbor means ghost cell mirror means matching energies means zero gradient.
						var cellA = this.cells[edge.cells[0]];
						var cellB = this.cells[edge.cells[1]];
						var dPressure = cellA.pressure - cellB.pressure;	//direction of normal
						var dx = cellA.x[0] - cellB.x[0];
						var dy = cellA.x[1] - cellB.x[1];
						var ds = Math.sqrt(dx * dx + dy * dy);
						
						cellA.q[1] -= dt * dPressure / ds * edge.normal[0];
						cellA.q[2] -= dt * dPressure / ds * edge.normal[1];
						
						cellB.q[1] += dt * dPressure / ds * edge.normal[0];
						cellB.q[2] += dt * dPressure / ds * edge.normal[1];
					}
				}
				/**/
				/*for (var i = 0; i < this.cells.length; ++i) {
					var cell = this.cells[i];
					var total_dPressure_times_dx = 0;
					var total_dPressure_times_dy = 0;
					var total_dx = 0;
					var total_dy = 0;
					for (var j = 0; j < cell.edges.length; ++j) {
						var edge = this.edges[cell.edges[j]];
						if (edge.cells.length != 2) continue;
						var cellB = edge.getOpposing(cell);
						var this_dPressure = cell.pressure - cellB.pressure;
						var dx = cell.x[0] - cellB.x[0];
						var dy = cell.x[1] - cellB.x[1];
						
						if (dx < 0) {
							total_dx -= dx;
							total_dPressure_times_dx -= this_dPressure;
						} else {
							total_dx += dx;
							total_dPressure_times_dx += this_dPressure;
						}
						
						if (dy < 0) {
							total_dy -= dy;
							total_dPressure_times_dy -= this_dPressure;
						} else {
							total_dy += dy;
							total_dPressure_times_dy += this_dPressure;
						}
					}
					var total_dPressure_dx = total_dPressure_times_dx / total_dx;
					var total_dPressure_dy = total_dPressure_times_dy / total_dy;
					
					cell.q[1] -= dt * total_dPressure_dx;
					cell.q[2] -= dt * total_dPressure_dy;
				}
				*/

				//apply work diffusion = momentum
				for (var i = 0; i < this.edges.length; ++i) {
					var edge = this.edges[i];
					if (edge.cells.length == 1) {
						var cellA = this.cells[edge.cells[0]];

						var dx = 2 * (cellA.x[0] - edge.x[0]);
						var dy = 2 * (cellA.x[1] - edge.x[1]);
						var ds = Math.sqrt(dx * dx + dy * dy);

						var rhoA = cellA.q[0];
						var uA = cellA.q[1] / rhoA;
						var vA = cellA.q[2] / rhoA;
						var pressureA = cellA.pressure;
						
						var rhoB = rhoA;
						var uB = -uA;
						var vB = -vA;
						var pressureB = pressureA;	
						var dWork = pressureA * (uA * edge.normal[0] + vA * edge.normal[1])
								- pressureB * (uB * edge.normal[0] + vB * edge.normal[1]);

						cellA.q[3] -= dt * dWork / ds; 
					} else {
						var cellA = this.cells[edge.cells[0]];
						var cellB = this.cells[edge.cells[1]];
						
						var dx = cellA.x[0] - cellB.x[0];
						var dy = cellA.x[1] - cellB.x[1];
						var ds = Math.sqrt(dx * dx + dy * dy);
				
						var rhoA = cellA.q[0];
						var uA = cellA.q[1] / rhoA;
						var vA = cellA.q[2] / rhoA;
						
						var rhoB = cellB.q[0];
						var uB = cellB.q[1] / rhoB;
						var vB = cellB.q[2] / rhoB;

						var dWork = cellA.pressure * (uA * edge.normal[0] + vA * edge.normal[1])
								- cellB.pressure * (uB * edge.normal[0] + vB * edge.normal[1]);

						cellA.q[3] -= dt * dWork / ds; 
						cellB.q[3] += dt * dWork / ds; 
					}
				}

			},
			'Riemann / Roe' : function(dt) {
				var qLs = [];
				var qRs = [];
				for (var j = 1; j < this.nx; ++j) {
					for (var i = 1; i < this.nx; ++i) {
						for (var side = 0; side < 2; ++side) {
							var indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
							var indexR = i + this.nx * j;
							var interfaceIndex = side + 2 * (i + (this.nx+1) * j);
							var solidL = this.solid[indexL];
							var solidR = this.solid[indexR];
							for (var k = 0; k < 4; ++k) {
								var qL = this.q[k + 4 * indexL];
								var qR = this.q[k + 4 * indexR];
								if (solidL && !solidR) {
									qL = qR;
									if (k == side+1) qL = -qL;
								} else if (solidR && !solidL) {
									qR = qL;
									if (k == side+1) qR = -qR;
								}
								qLs[k] = qL;
								qRs[k] = qR;
							}
							for (var state = 0; state < 4; ++state) {
								//find state change across interface in the basis of the eigenspace at the interface
								var sum = 0;
								for (var k = 0; k < 4; ++k) {
										//reproject into interface eigenspace
									sum += this.interfaceEigenvectorsInverse[state + 4 * (k + 4 * interfaceIndex)]
										//flux difference
										* (qRs[k] - qLs[k]);
								}
								this.interfaceDeltaQTilde[state + 4 * interfaceIndex] = sum;
							}
						}
					}
				}
				
				//boundary zero
				for (var j = 0; j < this.nghost-1; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var state = 0; state < 4; ++state) {
							//left boundary, left and top sides, zero vector
							this.interfaceDeltaQTilde[state + 4 * (0 + 2 * (j + (this.nx+1) * i))] = 0;
							//right boundary, left and top sides, zero vector
							this.interfaceDeltaQTilde[state + 4 * (0 + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
							//top boundary, left and top sides, zero vector
							this.interfaceDeltaQTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * j))] = 0;
							//bottom boundary, left and top sides, zero vector
							this.interfaceDeltaQTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
						}
					}
				}

				for (var j = this.nghost; j < this.nx + this.nghost - 3; ++j) {
					for (var i = this.nghost; i < this.nx + this.nghost - 3; ++i) {
						for (var side = 0; side < 2; ++side) {
							var indexL2 = i - 2 * dirs[side][0] + this.nx * (j - 2 * dirs[side][1]);
							var indexL1 = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
							var indexR1 = i + this.nx * j;
							var indexR2 = i + dirs[side][0] + this.nx * (j + dirs[side][1]); 
							
							var solidL2 = this.solid[indexL2];
							var solidL1 = this.solid[indexL1];
							var solidR1 = this.solid[indexR1];
							var solidR2 = this.solid[indexR2];
							
							var interfaceIndexL = side + 2 * (i - dirs[side][0] + (this.nx+1) * (j - dirs[side][1]));
							var interfaceIndex = side + 2 * (i + (this.nx+1) * j);
							var interfaceIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
							
							for (var state = 0; state < 4; ++state) {
								
								var interfaceDeltaQTildeL = this.interfaceDeltaQTilde[state + 4 * interfaceIndexL];
								var interfaceDeltaQTilde = this.interfaceDeltaQTilde[state + 4 * interfaceIndex];
								var interfaceDeltaQTildeR = this.interfaceDeltaQTilde[state + 4 * interfaceIndexR];

								//TODO is this right? 
								if (solidL2) {
									interfaceDeltaQTildeL = interfaceDeltaQTilde;
								}
								if (solidR2) {
									interfaceDeltaQTildeR = interfaceDeltaQTilde;
								}

								if (Math.abs(interfaceDeltaQTilde) > 0) {
									if (this.interfaceEigenvalues[state + 4 * interfaceIndex] > 0) {
										this.rTilde[state + 4 * interfaceIndex] = 
											this.interfaceDeltaQTilde[state + 4 * interfaceIndexL]
											/ interfaceDeltaQTilde;
									} else {
										this.rTilde[state + 4 * interfaceIndex] = 
											this.interfaceDeltaQTilde[state + 4 * interfaceIndexR]
											/ interfaceDeltaQTilde;
									}
								} else {
									this.rTilde[state + 4 * interfaceIndex] = 0;
								}
							}
						}
					}
				}
		
				//..and keep the boundary rTilde's zero	
				for (var j = 0; j < this.nghost; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var state = 0; state < 4; ++state) {
							//left
							this.rTilde[state + 4 * (0 + 2 * (j + (this.nx+1) * i))] = 0;
							//right
							this.rTilde[state + 4 * (0 + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
							//bottom
							this.rTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * j))] = 0;
							//top
							this.rTilde[state + 4 * (1 + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
						}
					}
				}
			
				var fluxAvg = [];	//4
				var fluxTilde = [];	//4
				var dxi = [];
				//transform cell q's into cell qTilde's (eigenspace)
				for (var j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
					for (var i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						var dx = this.xi[0 + 2 * (i + (this.nx+1) * j)] 
							- this.xi[0 + 2 * (i-dirs[0][0] + (this.nx+1) * (j-dirs[0][1]))];
						var dy = this.xi[1 + 2 * (i + (this.nx+1) * j)] 
							- this.xi[1 + 2 * (i-dirs[1][0] + (this.nx+1) * (j-dirs[1][1]))];
						var volume = dx * dy;
						dxi[0] = dx;
						dxi[1] = dy;
						for (var side = 0; side < 2; ++side) {

							var indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
							var indexR = i + this.nx * j;
							var solidL = this.solid[indexL];
							var solidR = this.solid[indexR];
							var interfaceIndex = side + 2 * (i + (this.nx+1) * j);
							
							for (var k = 0; k < 4; ++k) {
								var qL = this.q[k + 4 * indexL];
								var qR = this.q[k + 4 * indexR];
								if (solidL && !solidR) {
									qL = qR;
									if (k == side+1) qL = -qL;
								} else if (solidR && !solidL) {
									qR = qL;
									if (k == side+1) qR = -qR;
								}
								qLs[k] = qL;
								qRs[k] = qR;
							}
							
							//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
							//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
							for (var state = 0; state < 4; ++state) {
								var sum = 0;
								for (var k = 0; k < 4; ++k) {
									sum += this.interfaceMatrix[state + 4 * (k + 4 * interfaceIndex)]
										* (qRs[k] + qLs[k]);
								}
								fluxAvg[state] = .5 * sum;
							}

							//calculate flux
							for (var state = 0; state < 4; ++state) {
								var theta = 0;
								var eigenvalue = this.interfaceEigenvalues[state + 4 * interfaceIndex];
								if (eigenvalue >= 0) {
									theta = 1;
								} else {
									theta = -1;
								}
							
								var phi = fluxMethods[this.fluxMethod](this.rTilde[state + 4 * interfaceIndex]);
								var epsilon = eigenvalue * dt / dxi[side];//* volume / (dxi[side] * dxi[side]); 
								var deltaFluxTilde = eigenvalue * this.interfaceDeltaQTilde[state + 4 * interfaceIndex];
								fluxTilde[state] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
							}
						
							//reproject fluxTilde back into q
							for (var state = 0; state < 4; ++state) {
								var sum = 0;
								for (var k = 0; k < 4; ++k) {
									sum += fluxTilde[k] * this.interfaceEigenvectors[state + 4 * (k + 4 * interfaceIndex)];
								}
								this.flux[state + 4 * interfaceIndex] = fluxAvg[state] + sum;
							}
						}
					}
				}
			
				//zero boundary flux
				//..and keep the boundary r's zero	
				for (var j = 0; j < this.nghost-1; ++j) {
					for (var i = 0; i <= this.nx; ++i) {
						for (var state = 0; state < 4; ++state) {
							//left
							this.flux[state + 4 * (0 + 2 * (j + (this.nx+1) * i))] = 0;
							//right
							this.flux[state + 4 * (0 + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
							//bottom
							this.flux[state + 4 * (1 + 2 * (i + (this.nx+1) * j))] = 0;
							//top
							this.flux[state + 4 * (1 + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
						}
					}
				}

				//update cells
				for (var j = this.nghost; j < this.nx-this.nghost; ++j) {
					for (var i = this.nghost; i < this.nx-this.nghost; ++i) {
						if (this.solid[i + this.nx * j]) continue;
						
						var xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
						var xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
						var dx = this.xi[xiIndexR] - this.xi[xiIndexL];
						
						var xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
						var xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
						var dy = this.xi[xiIndexR] - this.xi[xiIndexL];
						
						var volume = dx * dy;
						dxi[0] = dx;
						dxi[1] = dy;
						
						for (var side = 0; side < 2; ++side) {
							var interfaceIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
							var interfaceIndexL = side + 2 * (i + (this.nx+1) * j);
							for (var state = 0; state < 4; ++state) {
								var interfaceFluxIndexR = state + 4 * interfaceIndexR;
								var interfaceFluxIndexL = state + 4 * interfaceIndexL;
								var df = this.flux[interfaceFluxIndexR] - this.flux[interfaceFluxIndexL];
								this.q[state + 4 * (i + this.nx * j)] -= dt * df / dxi[side];//* volume / (dxi[side] * dxi[side]);
							}
						}
					}
				}
			}
		}
	}
};

//called with 'this' the HydroState
var advectMethods = {
	Burgers : {
		initStep : function() {},
		calcCFLTimestep : function() {
			var mindum = undefined;
			for (var i = 0; i < this.cells.length; ++i) {
				var cell = this.cells[i];
				var x = cell.x[0];
				var y = cell.x[1];
				var rho = cell.q[0];
				var u = cell.q[1] / rho;
				var v = cell.q[2] / rho; 
				var energyTotal = cell.q[3] / rho; 
				var energyKinetic = .5 * (u * u + v * v);
				var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
				var energyThermal = energyTotal - energyKinetic - energyPotential;
				var speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
				
				//TODO fixme ... make this match Roe and run it across interface velocities?
				var dx = (xmax - xmin) / this.nx;
				var dy = (ymax - ymin) / this.ny;
				
				var dum = dx / (speedOfSound + Math.abs(u));
				if (mindum === undefined || dum < mindum) mindum = dum;
				var dum = dy / (speedOfSound + Math.abs(v));
				if (mindum === undefined || dum < mindum) mindum = dum;
			}
			return this.cfl * mindum;
		}
	},
	'Riemann / Roe' : {
		initStep : function() {
			for (var j = 1; j < this.nx; ++j) {
				for (var i = 1; i < this.nx; ++i) {
					for (var side = 0; side < 2; ++side) {
						var xL = side == 0 ? this.xi[0 + 2 * (i + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
						var xR = side == 0 ? this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
						var yL = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * j)] : this.x[1 + 2 * (i + this.nx * j)];
						var yR = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] : this.x[1 + 2 * (i + this.nx * j)];
				
						var normalX = dirs[side][0];
						var normalY = dirs[side][1];

						var indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
						var qIndexL = 4 * indexL;
						var solidL = this.solid[indexL];
						var densityL = this.q[0 + qIndexL];
						var velocityXL = this.q[1 + qIndexL] / densityL;
						var velocityYL = this.q[2 + qIndexL] / densityL;
						var energyTotalL = this.q[3 + qIndexL] / densityL;
						
						var indexR = i + this.nx * j;
						var qIndexR = 4 * indexR;
						var solidR = this.solid[indexR];
						var densityR = this.q[0 + qIndexR];
						var velocityXR = this.q[1 + qIndexR] / densityR;
						var velocityYR = this.q[2 + qIndexR] / densityR;
						var energyTotalR = this.q[3 + qIndexR] / densityR;					
						
						if (solidL && !solidR) {	//right cell has a wall on the left
							densityL = densityR;
							//reflect along normal
							var velocityNR = normalX * velocityXR + normalY * velocityYR;
							velocityXL = velocityXR - 2 * normalX * velocityNR;
							velocityYL = velocityYR - 2 * normalY * velocityNR;
							energyTotalL = energyTotalR;
						} else if (!solidL && solidR) {	//left cell has a wall on the right
							densityR = densityL;
							//reflect along normal
							var velocityNL = normalX * velocityXL + normalY * velocityYL;
							velocityXR = velocityXL - 2 * normalX * velocityNL;
							velocityYR = velocityYL - 2 * normalY * velocityNL;
							energyTotalR = energyTotalL;
						}
					
						var energyKineticL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
						var energyPotentialL = (xL - xmin) * externalForceX + (yL - ymin) * externalForceY;
						var energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
						var pressureL = (this.gamma - 1) * densityL * energyThermalL;
						var speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
						var hTotalL = energyTotalL + pressureL / densityL;
						var roeWeightL = Math.sqrt(densityL);
					
						var energyKineticR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
						var energyPotentialR = (xR - xmin) * externalForceX + (yR - ymin) * externalForceY;
						var energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
						var pressureR = (this.gamma - 1) * densityR * energyThermalR;
						var speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
						var hTotalR = energyTotalR + pressureR / densityR;
						var roeWeightR = Math.sqrt(densityR);

						var denom = roeWeightL + roeWeightR;
						var velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
						var velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
						var hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;
						
						buildEigenstate(
							 //index into interface element.  
							 //from there you'll have to scale by cell size.  
							 //Thus manually recreating the automatic storage of C structures. 
							 //JavaScript, why can't you be more like LuaJIT? 
							 side + 2 * (i + (this.nx+1) * j),	
							 this.interfaceMatrix,	//dim^2 = 16
							 this.interfaceEigenvalues,	//dim = 4
							 this.interfaceEigenvectors,	//dim^2 = 16
							 this.interfaceEigenvectorsInverse,	//dim^2 = 16
							 velocityX, velocityY, hTotal, this.gamma,
							 dirs[side][0], dirs[side][1]);
					}
				}
			}
		},
		calcCFLTimestep : function() {
			var mindum = undefined;
			for (var j = 1; j < this.nx; ++j) {
				for (var i = 1; i < this.nx; ++i) {
					if (this.solid[i + this.nx * j]) continue;
					for (var side = 0; side < 2; ++side) {
						var maxLambda = Math.max(0, 
							this.interfaceEigenvalues[0+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[1+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[2+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[3+4*(side+2*(i+(this.nx+1)*j))]);
						var minLambda = Math.min(0, 
							this.interfaceEigenvalues[0+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
							this.interfaceEigenvalues[1+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
							this.interfaceEigenvalues[2+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
							this.interfaceEigenvalues[3+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))]);
						var dx = this.xi[side + 2 * (i+dirs[side][0] + (this.nx+1) * (j+dirs[side][1]))] 
							- this.xi[side + 2 * (i + (this.nx+1) * j)];
						var dum = dx / (maxLambda - minLambda);
						if (mindum === undefined || dum < mindum) mindum = dum;
					}
				}
			}
			return this.cfl * mindum;
		}
	}
};

var HydroCell = function(args) {
	this.state = args.state;
	this.edges = args.edges;
	this.q = new Float64Array(4);	//state: rho, rho u, rho v, rho e
	this.pressure = 0;

	this.x = [0,0];
	for (var i = 0; i < this.edges.length; ++i) {
		var edge = this.state.edges[this.edges[i]];
		for (var j = 0; j < 2; ++j) {
			var vertex = this.state.vertexes[edge.vertexes[j]];
			this.x[0] += vertex[0];
			this.x[1] += vertex[1];
		}
	}
	this.x[0] /= 2 * this.edges.length;
	this.x[1] /= 2 * this.edges.length;
};

var HydroEdge = function(args) {
	this.state = args.state;
	this.vertexes = args.vertexes;
	
	var va = this.state.vertexes[this.vertexes[0]];
	var vb = this.state.vertexes[this.vertexes[1]];

	this.x = [
		.5 * (va[0] + vb[0]),
		.5 * (va[1] + vb[1])
	];

	var dx = vb[0] - va[0];
	var dy = vb[1] - va[1];
	//normalize
	var ds = Math.sqrt(dx * dx + dy * dy);
	dx /= ds;
	dy /= ds;
	
	this.normal = [-dy, dx];
	
	this.cells = [];


	//flux per-state
	this.flux = new Float64Array(4);


	//used for Burgers:


	//interface velocity
	this.u = 0;

	//used with slope limiters
	this.r = new Float64Array(4);


	//used for Riemann
/*

	this.jacobian = new Float64Array(4 * 4);
	this.eigenvalues = new Float64Array(4);
	this.eigenvectors = new Float64Array(4 * 4);
	this.eigenvectorsInverse = new Float64Array(4 * 4);
	for (var j = 0; j < 4; ++j) {
		this.eigenvalues[i] = 1;
		for (var i = 0; i < 4; ++i) {
			var delta_ij = i == j ? 1 : 0;
			this.jacobian[i + 4 * j] = delta_ij;
			this.eigenvectors[i + 4 * j] = delta_ij;
			this.eigenvetorsInverse[i + 4 * j] = delta_ij
		}
	}

	this.rTilde = new Float64Array(4);
*/

}
HydroEdge.prototype = {
	//after cells have added to boundaries, line up cells[0] with normal 
	initAux : function() {

		//also rearrange cells so that the normal points towards cells[0] and away from cells[1]
		if (this.cells.length < 1) throw "too few cells";
		if (this.cells.length > 2) throw "too many cells";
		if (this.cells.length == 2) {
			var cellAX = this.state.cells[this.cells[0]].x;
			var cellBX = this.state.cells[this.cells[1]].x;
			var dax = [cellAX[0] - this.x[0], cellAX[1] - this.x[1]];
			var dbx = [cellBX[0] - this.x[0], cellBX[1] - this.x[1]];
			var aforward = dax[0] * this.normal[0] + dax[1] * this.normal[1] > 0;
			var bforward = dbx[0] * this.normal[0] + dbx[1] * this.normal[1] > 0;
			if (aforward && bforward) throw "both cells in front of normal";
			if (!aforward && !bforward) throw "both cells behind normal";
			if (!aforward && bforward) {
				var tmp = this.cells[0];
				this.cells[0] = this.cells[1];
				this.cells[1] = tmp;
			}
		}	
	},
	getOpposing(cell) {
		if (this.state.cells[this.cells[0]] == cell) return this.state.cells[this.cells[1]];
		if (this.state.cells[this.cells[1]] == cell) return this.state.cells[this.cells[0]];
	}
};

var HydroState = makeClass({ 
	init : function(args) {
		this.nx = args.size;
		this.cfl =.5;
		this.gamma = args.gamma;

		//define vertices
		
		this.vertexes = [];

		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				var x = i / this.nx * (xmax - xmin) + xmin;
				var y = j / this.nx * (ymax - ymin) + ymin;
				this.vertexes.push([x, y]);
			}
		}
		
		//define edges by vertices 
		
		this.edges = [];
		
		//x_{i-1/2},{j-1/2},dim: interface positions
		//0 <= i,j < this.nx+1
		//0 <= dim < 2
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					if ((side == 0 && i > 0) || (side == 1 && j > 0)) {
						var vertexes = [];
						vertexes.push(i + (this.nx+1) * j);
						vertexes.push(i - dirs[side][0] + (this.nx+1) * (j - dirs[side][1]));
						this.edges.push(new HydroEdge({
							state : this,
							vertexes : vertexes
						}));
					}
				}
			}
		}
		
		//define cells by their edges 
		
		this.cells = [];	

		//x_i,j,dim: cell positions
		//0 <= i < this.nx
		var thiz = this;
		for (var j = 0; j < this.nx; ++j) {
			for (var i = 0; i < this.nx; ++i) {
				var edges = [];
				var cellIndex = this.cells.length;
				var addEdgeWithVertexes = function(a,b) {
					for (var edgeIndex = 0; edgeIndex < thiz.edges.length; ++edgeIndex) {
						var edge = thiz.edges[edgeIndex];
						if ((edge.vertexes[0] == a && edge.vertexes[1] == b) ||
							(edge.vertexes[0] == b && edge.vertexes[1] == a)) 
						{
							edge.cells.push(cellIndex);
							edges.push(edgeIndex);
							return;
						}
					}
					throw "couldn't find that edge";
				};
				addEdgeWithVertexes(i + (this.nx + 1) * j, i + 1 + (this.nx + 1) * j);
				addEdgeWithVertexes(i + (this.nx + 1) * j, i + (this.nx + 1) * (j + 1));
				addEdgeWithVertexes(i + 1 + (this.nx + 1) * j, i + 1 + (this.nx + 1) * (j + 1));
				addEdgeWithVertexes(i + (this.nx + 1) * (j + 1), i + 1 + (this.nx + 1) * (j + 1));
				this.cells.push(new HydroCell({
					state : this,
					edges : edges
				}));
			}
		}

		//init after we've added cells to interfaces
		for (var i = 0; i < this.edges.length; ++i) {
			this.edges[i].initAux();
		}
	
		this.resetSod();
		
		//qiTilde_{i-1/2},{j-1/2},side,state	
		this.interfaceDeltaQTilde = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4);
		for (var j = 0; j <= this.nx; ++j) {
			for (var i = 0; i <= this.nx; ++i) {
				for (var side = 0; side < 2; ++side) {
					for (var state = 0; state < 4; ++state) {
						this.interfaceDeltaQTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
					}
				}
			}
		}

		//solver configuration
		this.fluxMethod = 'donor cell';
		this.advectMethod = 'Burgers';
		this.integrationMethod = 'Forward Euler';
		this.drawToScreenMethod = 'Density';
	},
	resetSod : function() {
		for (var i = 0; i < this.cells.length; ++i) {
			var cell = this.cells[i];
			var x = cell.x[0];
			var y = cell.x[1];
			var rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
			var u = 0;
			var v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			var energyKinetic = .5 * (u * u + v * v);
			var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			var energyThermal = 1;
			var energyTotal = energyKinetic + energyThermal + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u; 
			cell.q[2] = rho * v; 
			cell.q[3] = rho * energyTotal; 
			if (cell.q[0] !== cell.q[0]) throw 'nan';
		}
	},
	resetSodCylinder : function() {
		//TODO restructure the mesh to avoid the cylinder
		for (var i = 0; i < this.cells.length; ++i) {
			var cell = this.cells[i];
			var x = cell.x[0];
			var y = cell.x[1];
			var rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
			var u = 0;
			var v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			var energyKinetic = .5 * (u * u + v * v);
			var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			var energyThermal = 1;
			var energyTotal = energyKinetic + energyThermal + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u; 
			cell.q[2] = rho * v; 
			cell.q[3] = rho * energyTotal; 
		}
	},
	resetWave : function() {
		var xmid = .5 * (xmin + xmax);
		var ymid = .5 * (ymin + ymax);
		var dg = .2 * (xmax - xmin);
		for (var i = 0; i < this.cells.length; ++i) {
			var cell = this.cells[i];
			var x = cell.x[0];
			var y = cell.y[0];
			var dx = x - xmid;
			var dy = y - ymid;
			var rho = 3 * Math.exp(-(dx * dx + dy * dy) / (dg * dg)) + .1;
			var u = 0; 
			var v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			var energyKinetic = .5 * (u * u + v * v);
			var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			var energyThermal = 1;
			var energyTotal = energyKinetic + energyThermal + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u;
			cell.q[2] = rho * v;
			cell.q[3] = rho * energyTotal;
		}
	},
	//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
	resetKelvinHemholtz : function() {
		var xmid = .5 * (xmin + xmax);
		for (var i = 0; i < this.cells.length; ++i) {
			var cell = this.cells[i];
			var x = cell.x[0];
			var y = cell.x[1];
			var yInTheMiddle = y > (.75 * ymin + .25 * ymax) && y < (.25 * ymin + .75 * ymax);
			var rho = yInTheMiddle ? 2 : 1;
			var u = yInTheMiddle ? .5 : -.5;
			var v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			//P = (gamma - 1) rho (eTotal - eKinetic - ePotential)
			//eTotal = P / ((gamma - 1) rho) + eKinetic + ePotential
			//eTotal = eThermal + eKinetic + ePotential
			//eThermal = P / ((gamma - 1) rho)
			var pressure = 2.5;
			var energyKinetic = .5 * (u * u + v * v);
			var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			var energyTotal = pressure / ((this.gamma - 1) * rho) + energyKinetic + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u; 
			cell.q[2] = rho * v;
			cell.q[3] = rho * energyTotal; 
		}
		//TODO make it periodic on the left/right borders and reflecting on the top/bottom borders
	},
	//http://www.astro.virginia.edu/VITA/ATHENA/rt.html
	resetRayleighTaylor : function() {
		var ymid = .5 * (ymin + ymax);
		for (var i = 0; i < this.cells.length; ++i) {
			var cell = this.cells[i];
			var x = cell.x[0 + 2 * e];
			var y = cell.x[1 + 2 * e];
			var yGreaterThanMid = y > ymid;
			var rho = yGreaterThanMid ? 2 : 1;
			var u = 0;
			var v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			//ePotential = g y
			//eThermal = eTotal - eKinetic - ePotential
			//P = (gamma - 1) rho eThermal = (gamma - 1) rho (eTotal - eKinetic - ePotential)
			//eTotal = P / ((gamma - 1) rho) + eKinetic + ePotential
			//...but Bernoulli's equation says
			//constant = P / rho * gamma / (gamma - 1) + eKinetic + ePotential
			var energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			var pressure = 2.5 - rho * energyPotential; 
			//var pressure = (this.gamma - 1) * rho * (2.5 - energyPotential);	//this fits with our non-external-force static ...
			/*
			pressure gradient zero ...
			dPressure / dy + rho * externalForce[side] = 0
			...integrate wrt y ...
			pressure + rho * externalForce.y * y = constant
			pressure = constant - rho * externalForce.y = constant - rho * energyPotential
			*/
			var energyKinetic = .5 * (u * u + v * v);
			/*
			energyTotal = energyKinetic + energyPotential + energyThermal
			energyTotal = .5 * (u*u + v*v) + energyPotential + pressure / (rho * (gamma - 1))
			*/
			var energyTotal = pressure / (rho * (this.gamma - 1)) + energyKinetic + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u;
			cell.q[2] = rho * v;
			cell.q[3] = rho * energyTotal;
		}
	},
	step : function(dt) {
		//solve
		integrationMethods[this.integrationMethod].advect[this.advectMethod].call(this, dt);
	},
	update : function() {
//console.log('update');
		//do any pre-calcCFLTimestep preparation (Roe computes eigenvalues here)
		advectMethods[this.advectMethod].initStep.call(this);
		
		//get timestep
		var dt;
		if (useCFL) {
			dt = advectMethods[this.advectMethod].calcCFLTimestep.call(this);
		} else {
			dt = fixedDT;
		}
window.lastDT = dt;
//if (dt !== dt) throw 'nan';
		//do the update
		this.step(dt);
	}
});


var Hydro = makeClass({
	init : function() {
		var size = Number($.url().param('size'));
		if (size === undefined || size !== size) size = 50;
		var gamma = Number($.url().param('gamma'));
		if (gamma === undefined || gamma !== gamma) gamma = 7/5;
		this.state = new HydroState({
			size : size,
			gamma : gamma
		});

		var plainShader = new GL.ShaderProgram({
			vertexCode : mlstr(function(){/*
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
	gl_PointSize = 2.;
}
*/}),
			vertexPrecision : 'best',
			fragmentCode : mlstr(function(){/*
void main() {
	gl_FragColor = vec4(1.);
}
*/}),
			fragmentPrecision : 'best'
		});

		var stateShader = new GL.ShaderProgram({
			vertexCode : mlstr(function(){/*
attribute vec2 vertex;
attribute float state;
varying float statev;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	statev = state;
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
	gl_PointSize = 3.;
}
*/}),
			vertexPrecision : 'best',
			fragmentCode : mlstr(function(){/*
varying float statev;
uniform sampler2D tex;
void main() {
	gl_FragColor = texture2D(tex, vec2(statev, .5));
}
*/}),
			fragmentPrecision : 'best'
		});

		//vertices
		var vertexData = [];
		for (var i = 0; i < this.state.vertexes.length; ++i) {
			for (var j = 0; j < 2; ++j) {
				vertexData.push(this.state.vertexes[i][j]);
			}
		}
		var vertexBuffer = new GL.ArrayBuffer({dim : 2, data : vertexData, keep : true});
		this.vertexObj = new GL.SceneObject({
			mode : gl.POINTS,
			attrs : {
				vertex : vertexBuffer
			},
			shader : plainShader,
			static : true
		});
	
		//edges
		var indexData = [];
		for (var i = 0; i < this.state.edges.length; ++i) {
			for (var j = 0; j < 2; ++j) {
				indexData.push(this.state.edges[i].vertexes[j]);
			}
		}
		var indexBuffer = new GL.ElementArrayBuffer({data : indexData})
		this.edgeObj = new GL.SceneObject({
			mode : gl.LINES,
			indexes : indexBuffer,
			attrs : {
				vertex : vertexBuffer
			},
			shader : plainShader,
			static : true
		});
	
		//cells
		//points for now
		//I should make polys out of the cells or something
		var cellVertexData = [];
		var cellStateData = [];
		for (var i = 0; i < this.state.cells.length; ++i) {
			var cell = this.state.cells[i];
			for (var j = 0; j < 2; ++j) {
				cellVertexData.push(cell.x[j]);
			}
			cellStateData.push(0);
		}
		var cellVertexBuffer = new GL.ArrayBuffer({dim : 2, data : cellVertexData, keep : true});
		var cellStateBuffer = new GL.ArrayBuffer({dim : 1, data : cellStateData, keep : true});
		this.cellObj = new GL.SceneObject({
			mode : gl.POINTS,
			attrs : {
				vertex : cellVertexBuffer,
				state : cellStateBuffer
			},
			shader : stateShader,
			uniforms : {
				tex : 0
			},
			static : true
		});
	},
	update : function() {
		//todo adm or something
		//update a copy of the grid and its once-refined
		//...and a once-unrefined ... over mergeable cells only?
		//then test for errors and split when needed
		this.state.update();

		//update geometry z coordinate
		var result = drawToScreenMethods[this.state.drawToScreenMethod].call(this);
		var dataMin = result[0];
		var dataMax = result[1];
		
		if (this.updateLastDataRange) {
			this.lastDataMin = dataMin;
			this.lastDataMax = dataMax;
			var thisTick = Math.floor(Date.now() / 1000);
			if (this.lastDataRangeUpdateTime != thisTick) {
				this.lastDataRangeUpdateTime = thisTick;
				$('#dataRangeFixedMin').val(this.lastDataMin);
				$('#dataRangeFixedMax').val(this.lastDataMax);
			}
		}
		
		//update state buffers
		this.cellObj.attrs.state.updateData();
	}
});

var hydro;
function update() {
	//iterate
	hydro.update();
	//draw
	GL.draw();
	requestAnimFrame(update);
}

function onresize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	GL.resize();
}

function buildSelect(id, key, map) {
	var select = $('#' + id);
	for (var k in map) {
		var option = $('<option>', {text : k});
		option.appendTo(select);
		if (hydro.state[key] == k) {
			option.attr('selected', 'true');
		}
	}
	select.change(function() {
		hydro.state[key] = select.val();
	});
	return select;
}

var colorSchemes = {};

$(document).ready(function(){
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()
	
	try {
		gl = GL.init(canvas);
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}

	hydro = new Hydro();


	
	panel = $('#panel');	

	var panelContent = $('#content');
	$('#menu').click(function() {
		if (panelContent.css('display') == 'none') {
			panelContent.show();
		} else {
			panelContent.hide();
		}
	});

	$('#resetSod').click(function(){ hydro.state.resetSod(); });
	$('#resetSodCylinder').click(function(){ hydro.state.resetSodCylinder(); });
	$('#resetWave').click(function(){ hydro.state.resetWave(); });
	$('#resetKelvinHemholtz').click(function(){ hydro.state.resetKelvinHemholtz(); });
	$('#resetRayleighTaylor').click(function(){ hydro.state.resetRayleighTaylor(); });

	$('#useNoise').change(function() {
		useNoise = $(this).is(':checked');
	});

	(function(){
		var select = $('#gridsize');
		$.each([20, 50, 100, 200, 500, 1000], function(i,gridsize){
			var option = $('<option>', {text : gridsize});
			if (hydro.state.nx == gridsize) option.attr('selected', 'true');
			option.appendTo(select);
		});
		select.change(function(){
			var params = $.url().param();
			params.size = select.val();
			var url = location.href.match('[^?]*');
			var sep = '?';
			for (k in params) {
				if (k == '' && params[k] == '') continue;
				url += sep;
				url += k + '=' + params[k];
				sep = '&';
			}
			location.href = url;
		});
	})();

	buildSelect('fluxMethod', 'fluxMethod', fluxMethods);
	buildSelect('advectMethod', 'advectMethod', advectMethods);
	buildSelect('integrationMethod', 'integrationMethod', integrationMethods);
	buildSelect('drawToScreenMethod', 'drawToScreenMethod', drawToScreenMethods);

	$.each([
		'externalForceX',
		'externalForceY'
	], function(i,varName) {
		$('#'+varName).val(window[varName]);
		$('#'+varName).change(function() {
			var v = Number($(this).val());
			if (v !== v) return;	//NaN
			window[varName] = v;
		});
	});
	
	hydro.lastDataMin = Number($('#dataRangeFixedMin').val());
	hydro.lastDataMax = Number($('#dataRangeFixedMax').val());
	hydro.updateLastDataRange = true;
	$('#dataRangeScaleNormalized').change(function() {
		if (!$(this).is(':checked')) return;
		hydro.updateLastDataRange = true;
	});
	$('#dataRangeScaleFixed').change(function() {
		if (!$(this).is(':checked')) return;
		hydro.updateLastDataRange = false;
		hydro.lastDataMin = Number($('#dataRangeFixedMin').val()); 
		hydro.lastDataMax = Number($('#dataRangeFixedMax').val()); 
	});
	$('#dataRangeFixedMin').change(function() {
		if (hydro.updateLastDataRange) return;
		hydro.lastDataMin = Number($('#dataRangeFixedMin').val()); 
	});
	$('#dataRangeFixedMax').change(function() {
		if (hydro.updateLastDataRange) return;
		hydro.lastDataMax = Number($('#dataRangeFixedMax').val()); 
	});

	$('#timeStepCFLBased').change(function() {
		if (!$(this).is(':checked')) return;
		useCFL = true;
	});
	$('#timeStepCFL').val(hydro.state.cfl);
	$('#timeStepCFL').change(function() {
		var v = Number($(this).val());
		if (v != v) return;
		hydro.state.cfl = v;
	});
	$('#timeStepFixed').change(function() {
		if (!$(this).is(':checked')) return;
		useCFL = false;
	});
	$('#timeStepValue').val(fixedDT);
	$('#timeStepValue').change(function() {
		var v = Number($(this).val());
		if (v != v) return;	//stupid javascript ... convert anything it doesn't understand to NaNs...
		fixedDT = v;
	});

	gl.enable(gl.DEPTH_TEST);
	GL.view.ortho = true;
	GL.view.zNear = -1;
	GL.view.zFar = 1;
	GL.view.fovY = 125 / 200 * (xmax - xmin);

	colorSchemes.Heat = new GL.GradientTexture({
		width:256, 
		colors:[
			[0,0,0],
			[0,0,1],
			[1,1,0],
			[1,0,0],
		],
		dontRepeat : true
	});

	var isobarSize = 16;
	var isobarData = new Uint8Array(isobarSize);
	for (var i = 1; i < isobarSize; i += 2) {
		isobarData[i] = 255;
	}
	colorSchemes['B&W'] = new GL.Texture2D({
		width : isobarSize,
		height : 1,
		format : gl.LUMINANCE,
		internalFormat : gl.LUMINANCE,
		data : isobarData,
		minFilter : gl.LINEAR,
		magFilter : gl.NEAREST,
		generateMipmap : true,
		wrap : {s : gl.REPEAT, t : gl.CLAMP_TO_EDGE }
	});

	colorSchemes.Heat.bind();

	for (var colorSchemeName in colorSchemes) {
		(function(){
			var v = colorSchemes[colorSchemeName];
			$('#colorScheme').append($('<option>', {
				value : colorSchemeName,
				text : colorSchemeName
			}));
		})();
	}
	$('#colorScheme').change(function() {
		var colorSchemeName = $(this).val();
		gl.bindTexture(gl.TEXTURE_2D, colorSchemes[colorSchemeName].obj);
	});

	var zoomFactor = .0003;	// upon mousewheel
	var dragging = false;
	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
			dragging = false;
		},
		move : function(dx,dy) {
			dragging = true;
			var aspectRatio = canvas.width / canvas.height;
			GL.view.pos[0] -= dx / canvas.width * 2 * (aspectRatio * GL.view.fovY);
			GL.view.pos[1] += dy / canvas.height * 2 * GL.view.fovY;
			GL.updateProjection();
		},
		zoom : function(zoomChange) {
			dragging = true;
			var scale = Math.exp(-zoomFactor * zoomChange);
			GL.view.fovY *= scale 
			GL.updateProjection();
		}
	});
	
	//start it off
	$(window).resize(onresize);
	onresize();
	update();
});
