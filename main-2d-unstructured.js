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

let panel;
let canvas;

let xmin = -.5;
let xmax = .5; 
let ymin = -.5;
let ymax = .5;

let useNoise = true;
let noiseAmplitude = .01;

let useCFL = true;
let fixedDT = .2;
let gaussSeidelIterations = 20;

let mouse;

let externalForceX = 0;
let externalForceY = 0;

//interface directions
let dirs = [[1,0], [0,1]];

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
	let speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)));
	let tangentX = -normalY;
	let tangentY = normalX;
	let velocityN = velocityX * normalX + velocityY * normalY;
	let velocityT = velocityX * tangentX + velocityY * tangentY;
	let velocitySq = velocityX * velocityX + velocityY * velocityY;	
	
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
	let identCheck = [];
	let identBad = false;
	for (let i = 0; i < 4; ++i) {
		for (let j = 0; j < 4; ++j) {
			let s = 0;
			let d = 0;
			for (let k = 0; k < 4; ++k) {
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
			let epsilon = 1e-5;
			if (Math.abs(d - (i == j ? 1 : 0)) > epsilon) identBad = true;
		}
	}
	if (identBad) {
		console.log('bad eigen basis', identCheck);
	}
	/** /	
	function f32subset(a, o, s) {
		let d = new Float64Array(s);
		for (let i = 0; i < s; ++i) {
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

let drawToScreenMethods = {
	Density : function() {
		let dataMin = this.state.cells[0].q[0];
		let dataMax = dataMin + 1e-9;
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		for (let i = 0; i < this.state.cells.length; ++i) {	
			let s = this.state.cells[i].q[0];
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Velocity : function() {
		let rho0 = this.state.cells[0].q[0];
		let u0 = this.state.cells[0].q[1] / rho0;
		let v0 = this.state.cells[0].q[2] / rho0;
		let dataMin = Math.sqrt(u0 * u0 + v0 * v0);
		let dataMax = dataMin + 1e-9; 
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		for (let i = 0; i < this.state.cells.length; ++i) {
			let cell = this.state.cells[i];
			let rho = cell.q[0];
			let mx = cell.q[1];
			let my = cell.q[2];
			let s = Math.sqrt(mx * mx + my * my) / rho;
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Energy : function() {
		let q0 = this.state.cells[0].q;
		let dataMin = q0[3] / q0[0]; 
		let dataMax = dataMin + 1e-9; 
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		for (let i = 0; i < this.state.cells.length; ++i) {
			let cell = this.state.cells[i];
			let s = cell.q[3] / cell.q[0];
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Pressure : function() {
		let pressure = this.state.cells[0].pressure;
		let dataMin = pressure[0];
		let dataMax = dataMin + 1e-9;
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		for (let i = 0; i < this.state.cells.length; ++i) {
			let s = this.state.cells[i].pressure;
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.cellObj.attrs.state.data[i] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	}/*,
	Curl : function() {
		let nx = this.state.nx;
		let q = this.state.q;
		let dataMin = undefined; 
		let dataMax = undefined; 
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		let e = 0;
		let dx = (xmax - xmin) / nx;
		let dy = (ymax - ymin) / nx;
		for (let j = 0; j < nx; ++j) {
			for (let i = 0; i < nx; ++i) {
				let rho = q[0 + 4 * e];
				let mx = q[1 + 4 * e];
				let my = q[2 + 4 * e];
				let mxPrev = this.state.q[1 + 4 * ( ((i+nx-1)%nx) + nx * j )];
				let myPrev = this.state.q[1 + 4 * ( i + nx * ((j+nx-1)%nx) )];
				let dmx = mx - mxPrev;
				let dmy = my - myPrev;
				let s = (dmx / dy - dmy / dx) / rho;
				if (dataMin === undefined || s < dataMin) dataMin = s;
				if (dataMax === undefined || s > dataMax) dataMax = s;
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		return [dataMin, dataMax];
	}*/
};

let integrationMethods = {
	'Forward Euler' : {
		advect : {
			Burgers : function(dt) {
				/**/
				for (let i = 0; i < this.edges.length; ++i) {
					//TODO weight averages by distances from edge ... or volume ... or something
					let edge = this.edges[i];
					if (edge.cells.length == 1) {
						//zero boundaries
						let cell = this.cells[edge.cells[0]];
						edge.u = 0;
					} else {
						let cellA = this.cells[edge.cells[0]];
						let cellB = this.cells[edge.cells[1]];
						let rhoA = cellA.q[0];
						let uA = cellA.q[1] / rhoA;
						let vA = cellA.q[2] / rhoA;
						let rhoB = cellB.q[0];
						let uB = cellB.q[1] / rhoB;
						let vB = cellB.q[2] / rhoB;
						edge.u = .5 * ((uA + uB) * edge.normal[0] + (vA + vB) * edge.normal[1]);
					}
					//if (edge.u !== edge.u) throw 'nan';
				}

				//compute flux and advect for each state vector
				for (let i = 0; i < this.edges.length; ++i) {
					let edge = this.edges[i];
					//TODO slope calc
					for (let j = 0; j < 4; ++j) {
						edge.r[j] = 0;
					}
				}

				//construct flux:
				for (let i = 0; i < this.edges.length; ++i) {
					let edge = this.edges[i];
					if (edge.cells.length == 1) {
						//TODO boundary reflection?
						edge.flux[0] = edge.flux[1] = edge.flux[2] = edge.flux[3] = 0;
					} else {
						let cellA = this.cells[edge.cells[0]];
						let cellB = this.cells[edge.cells[1]];
						if (edge.u >= 0) {	//pull from cell B
							for (let j = 0; j < 4; ++j) {
								edge.flux[j] = edge.u * cellB.q[j];
							}
						} else {
							for (let j = 0; j < 4; ++j) {
								edge.flux[j] = edge.u * cellA.q[j];
							}
						}
						//TODO add flux limiter stuff
					}
				}

				//if (dt !== dt) throw 'nan';
				for (let i = 0; i < this.cells.length; ++i) {
					let cell = this.cells[i];
					for (let j = 0; j < cell.edges.length; ++j) {
						let edge = this.edges[cell.edges[j]];
						let dx = edge.x[0] - cell.x[0];
						let dy = edge.x[1] - cell.x[1];
						let ds = Math.sqrt(dx * dx + dy * dy);
						//if (ds !== ds) throw 'nan';
						if (ds === 0) throw 'divide by zero for cell ' + i + ' and edge ' + cell.edges[j];
						for (let k = 0; k < 4; ++k) {
							//if (cell.q[k] !== cell.q[k]) throw 'nan';
							//if (edge.flux[k] !== edge.flux[k]) throw 'nan';
							cell.q[k] -= dt / ds * edge.flux[k];
							//if (cell.q[k] !== cell.q[k]) throw 'nan';
						}
					}
				}
				/**/
				
				//compute pressure
				for (let i = 0; i < this.cells.length; ++i) {
					let cell = this.cells[i];
					//if (cell.q[0] !== cell.q[0]) throw 'nan';
					//if (cell.q[1] !== cell.q[1]) throw 'nan';
					//if (cell.q[2] !== cell.q[2]) throw 'nan';
					//if (cell.q[0] === 0) throw 'divide by zero error';
					let x = cell.x[0];
					let y = cell.x[1];
					let rho = cell.q[0];
					let u = cell.q[1] / rho;
					let v = cell.q[2] / rho;
					let energyTotal = cell.q[3] / rho;
					let energyKinetic = .5 * (u * u + v * v);
					let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
					let energyThermal = energyTotal - energyKinetic - energyPotential;
					cell.pressure = (this.gamma - 1) * rho * energyThermal;
					//if (cell.pressure !== cell.pressure) throw 'nan';
				}

				//apply external force
				for (let i = 0; i < this.cells.length; ++i) {
					let cell = this.cells[i];
					let rho = cell.q[0];
					cell.q[3] -= dt * (externalForceX * cell.q[1] + externalForceY * cell.q[2]);
					cell.q[1] -= dt * rho * externalForceX;
					cell.q[2] -= dt * rho * externalForceY;
				}

				//apply momentum diffusion = pressure
				/* per-edge method */
				for (let i = 0; i < this.edges.length; ++i) {
					let edge = this.edges[i];
					if (edge.cells.length == 2) {	//one neighbor means ghost cell mirror means matching energies means zero gradient.
						let cellA = this.cells[edge.cells[0]];
						let cellB = this.cells[edge.cells[1]];
						let dPressure = cellA.pressure - cellB.pressure;	//direction of normal
						let dx = cellA.x[0] - cellB.x[0];
						let dy = cellA.x[1] - cellB.x[1];
						let ds = Math.sqrt(dx * dx + dy * dy);
						
						cellA.q[1] -= dt * dPressure / ds * edge.normal[0];
						cellA.q[2] -= dt * dPressure / ds * edge.normal[1];
						
						cellB.q[1] += dt * dPressure / ds * edge.normal[0];
						cellB.q[2] += dt * dPressure / ds * edge.normal[1];
					}
				}
				/**/
				/*for (let i = 0; i < this.cells.length; ++i) {
					let cell = this.cells[i];
					let total_dPressure_times_dx = 0;
					let total_dPressure_times_dy = 0;
					let total_dx = 0;
					let total_dy = 0;
					for (let j = 0; j < cell.edges.length; ++j) {
						let edge = this.edges[cell.edges[j]];
						if (edge.cells.length != 2) continue;
						let cellB = edge.getOpposing(cell);
						let this_dPressure = cell.pressure - cellB.pressure;
						let dx = cell.x[0] - cellB.x[0];
						let dy = cell.x[1] - cellB.x[1];
						
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
					let total_dPressure_dx = total_dPressure_times_dx / total_dx;
					let total_dPressure_dy = total_dPressure_times_dy / total_dy;
					
					cell.q[1] -= dt * total_dPressure_dx;
					cell.q[2] -= dt * total_dPressure_dy;
				}
				*/

				//apply work diffusion = momentum
				for (let i = 0; i < this.edges.length; ++i) {
					let edge = this.edges[i];
					if (edge.cells.length == 1) {
						let cellA = this.cells[edge.cells[0]];

						let dx = 2 * (cellA.x[0] - edge.x[0]);
						let dy = 2 * (cellA.x[1] - edge.x[1]);
						let ds = Math.sqrt(dx * dx + dy * dy);

						let rhoA = cellA.q[0];
						let uA = cellA.q[1] / rhoA;
						let vA = cellA.q[2] / rhoA;
						let pressureA = cellA.pressure;
						
						let rhoB = rhoA;
						let uB = -uA;
						let vB = -vA;
						let pressureB = pressureA;	
						let dWork = pressureA * (uA * edge.normal[0] + vA * edge.normal[1])
								- pressureB * (uB * edge.normal[0] + vB * edge.normal[1]);

						cellA.q[3] -= dt * dWork / ds; 
					} else {
						let cellA = this.cells[edge.cells[0]];
						let cellB = this.cells[edge.cells[1]];
						
						let dx = cellA.x[0] - cellB.x[0];
						let dy = cellA.x[1] - cellB.x[1];
						let ds = Math.sqrt(dx * dx + dy * dy);
				
						let rhoA = cellA.q[0];
						let uA = cellA.q[1] / rhoA;
						let vA = cellA.q[2] / rhoA;
						
						let rhoB = cellB.q[0];
						let uB = cellB.q[1] / rhoB;
						let vB = cellB.q[2] / rhoB;

						let dWork = cellA.pressure * (uA * edge.normal[0] + vA * edge.normal[1])
								- cellB.pressure * (uB * edge.normal[0] + vB * edge.normal[1]);

						cellA.q[3] -= dt * dWork / ds; 
						cellB.q[3] += dt * dWork / ds; 
					}
				}

			},
			'Riemann / Roe' : function(dt) {
				let qLs = [];
				let qRs = [];
				for (let j = 1; j < this.nx; ++j) {
					for (let i = 1; i < this.nx; ++i) {
						for (let side = 0; side < 2; ++side) {
							let indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
							let indexR = i + this.nx * j;
							let interfaceIndex = side + 2 * (i + (this.nx+1) * j);
							let solidL = this.solid[indexL];
							let solidR = this.solid[indexR];
							for (let k = 0; k < 4; ++k) {
								let qL = this.q[k + 4 * indexL];
								let qR = this.q[k + 4 * indexR];
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
							for (let state = 0; state < 4; ++state) {
								//find state change across interface in the basis of the eigenspace at the interface
								let sum = 0;
								for (let k = 0; k < 4; ++k) {
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
				for (let j = 0; j < this.nghost-1; ++j) {
					for (let i = 0; i <= this.nx; ++i) {
						for (let state = 0; state < 4; ++state) {
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

				for (let j = this.nghost; j < this.nx + this.nghost - 3; ++j) {
					for (let i = this.nghost; i < this.nx + this.nghost - 3; ++i) {
						for (let side = 0; side < 2; ++side) {
							let indexL2 = i - 2 * dirs[side][0] + this.nx * (j - 2 * dirs[side][1]);
							let indexL1 = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
							let indexR1 = i + this.nx * j;
							let indexR2 = i + dirs[side][0] + this.nx * (j + dirs[side][1]); 
							
							let solidL2 = this.solid[indexL2];
							let solidL1 = this.solid[indexL1];
							let solidR1 = this.solid[indexR1];
							let solidR2 = this.solid[indexR2];
							
							let interfaceIndexL = side + 2 * (i - dirs[side][0] + (this.nx+1) * (j - dirs[side][1]));
							let interfaceIndex = side + 2 * (i + (this.nx+1) * j);
							let interfaceIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
							
							for (let state = 0; state < 4; ++state) {
								
								let interfaceDeltaQTildeL = this.interfaceDeltaQTilde[state + 4 * interfaceIndexL];
								let interfaceDeltaQTilde = this.interfaceDeltaQTilde[state + 4 * interfaceIndex];
								let interfaceDeltaQTildeR = this.interfaceDeltaQTilde[state + 4 * interfaceIndexR];

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
				for (let j = 0; j < this.nghost; ++j) {
					for (let i = 0; i <= this.nx; ++i) {
						for (let state = 0; state < 4; ++state) {
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
			
				let fluxAvg = [];	//4
				let fluxTilde = [];	//4
				let dxi = [];
				//transform cell q's into cell qTilde's (eigenspace)
				for (let j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
					for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						let dx = this.xi[0 + 2 * (i + (this.nx+1) * j)] 
							- this.xi[0 + 2 * (i-dirs[0][0] + (this.nx+1) * (j-dirs[0][1]))];
						let dy = this.xi[1 + 2 * (i + (this.nx+1) * j)] 
							- this.xi[1 + 2 * (i-dirs[1][0] + (this.nx+1) * (j-dirs[1][1]))];
						let volume = dx * dy;
						dxi[0] = dx;
						dxi[1] = dy;
						for (let side = 0; side < 2; ++side) {

							let indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
							let indexR = i + this.nx * j;
							let solidL = this.solid[indexL];
							let solidR = this.solid[indexR];
							let interfaceIndex = side + 2 * (i + (this.nx+1) * j);
							
							for (let k = 0; k < 4; ++k) {
								let qL = this.q[k + 4 * indexL];
								let qR = this.q[k + 4 * indexR];
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
							for (let state = 0; state < 4; ++state) {
								let sum = 0;
								for (let k = 0; k < 4; ++k) {
									sum += this.interfaceMatrix[state + 4 * (k + 4 * interfaceIndex)]
										* (qRs[k] + qLs[k]);
								}
								fluxAvg[state] = .5 * sum;
							}

							//calculate flux
							for (let state = 0; state < 4; ++state) {
								let theta = 0;
								let eigenvalue = this.interfaceEigenvalues[state + 4 * interfaceIndex];
								if (eigenvalue >= 0) {
									theta = 1;
								} else {
									theta = -1;
								}
							
								let phi = fluxMethods[this.fluxMethod](this.rTilde[state + 4 * interfaceIndex]);
								let epsilon = eigenvalue * dt / dxi[side];//* volume / (dxi[side] * dxi[side]); 
								let deltaFluxTilde = eigenvalue * this.interfaceDeltaQTilde[state + 4 * interfaceIndex];
								fluxTilde[state] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta));
							}
						
							//reproject fluxTilde back into q
							for (let state = 0; state < 4; ++state) {
								let sum = 0;
								for (let k = 0; k < 4; ++k) {
									sum += fluxTilde[k] * this.interfaceEigenvectors[state + 4 * (k + 4 * interfaceIndex)];
								}
								this.flux[state + 4 * interfaceIndex] = fluxAvg[state] + sum;
							}
						}
					}
				}
			
				//zero boundary flux
				//..and keep the boundary r's zero	
				for (let j = 0; j < this.nghost-1; ++j) {
					for (let i = 0; i <= this.nx; ++i) {
						for (let state = 0; state < 4; ++state) {
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
				for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
					for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
						if (this.solid[i + this.nx * j]) continue;
						
						let xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
						let xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
						let dx = this.xi[xiIndexR] - this.xi[xiIndexL];
						
						let xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
						let xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
						let dy = this.xi[xiIndexR] - this.xi[xiIndexL];
						
						let volume = dx * dy;
						dxi[0] = dx;
						dxi[1] = dy;
						
						for (let side = 0; side < 2; ++side) {
							let interfaceIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
							let interfaceIndexL = side + 2 * (i + (this.nx+1) * j);
							for (let state = 0; state < 4; ++state) {
								let interfaceFluxIndexR = state + 4 * interfaceIndexR;
								let interfaceFluxIndexL = state + 4 * interfaceIndexL;
								let df = this.flux[interfaceFluxIndexR] - this.flux[interfaceFluxIndexL];
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
let advectMethods = {
	Burgers : {
		initStep : function() {},
		calcCFLTimestep : function() {
			let mindum = undefined;
			for (let i = 0; i < this.cells.length; ++i) {
				let cell = this.cells[i];
				let x = cell.x[0];
				let y = cell.x[1];
				let rho = cell.q[0];
				let u = cell.q[1] / rho;
				let v = cell.q[2] / rho; 
				let energyTotal = cell.q[3] / rho; 
				let energyKinetic = .5 * (u * u + v * v);
				let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
				let energyThermal = energyTotal - energyKinetic - energyPotential;
				let speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
				
				//TODO fixme ... make this match Roe and run it across interface velocities?
				let dx = (xmax - xmin) / this.nx;
				let dy = (ymax - ymin) / this.ny;
				
				let dum = dx / (speedOfSound + Math.abs(u));
				if (mindum === undefined || dum < mindum) mindum = dum;
				let dum = dy / (speedOfSound + Math.abs(v));
				if (mindum === undefined || dum < mindum) mindum = dum;
			}
			return this.cfl * mindum;
		}
	},
	'Riemann / Roe' : {
		initStep : function() {
			for (let j = 1; j < this.nx; ++j) {
				for (let i = 1; i < this.nx; ++i) {
					for (let side = 0; side < 2; ++side) {
						let xL = side == 0 ? this.xi[0 + 2 * (i + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
						let xR = side == 0 ? this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
						let yL = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * j)] : this.x[1 + 2 * (i + this.nx * j)];
						let yR = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] : this.x[1 + 2 * (i + this.nx * j)];
				
						let normalX = dirs[side][0];
						let normalY = dirs[side][1];

						let indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
						let qIndexL = 4 * indexL;
						let solidL = this.solid[indexL];
						let densityL = this.q[0 + qIndexL];
						let velocityXL = this.q[1 + qIndexL] / densityL;
						let velocityYL = this.q[2 + qIndexL] / densityL;
						let energyTotalL = this.q[3 + qIndexL] / densityL;
						
						let indexR = i + this.nx * j;
						let qIndexR = 4 * indexR;
						let solidR = this.solid[indexR];
						let densityR = this.q[0 + qIndexR];
						let velocityXR = this.q[1 + qIndexR] / densityR;
						let velocityYR = this.q[2 + qIndexR] / densityR;
						let energyTotalR = this.q[3 + qIndexR] / densityR;					
						
						if (solidL && !solidR) {	//right cell has a wall on the left
							densityL = densityR;
							//reflect along normal
							let velocityNR = normalX * velocityXR + normalY * velocityYR;
							velocityXL = velocityXR - 2 * normalX * velocityNR;
							velocityYL = velocityYR - 2 * normalY * velocityNR;
							energyTotalL = energyTotalR;
						} else if (!solidL && solidR) {	//left cell has a wall on the right
							densityR = densityL;
							//reflect along normal
							let velocityNL = normalX * velocityXL + normalY * velocityYL;
							velocityXR = velocityXL - 2 * normalX * velocityNL;
							velocityYR = velocityYL - 2 * normalY * velocityNL;
							energyTotalR = energyTotalL;
						}
					
						let energyKineticL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
						let energyPotentialL = (xL - xmin) * externalForceX + (yL - ymin) * externalForceY;
						let energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
						let pressureL = (this.gamma - 1) * densityL * energyThermalL;
						let speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
						let hTotalL = energyTotalL + pressureL / densityL;
						let roeWeightL = Math.sqrt(densityL);
					
						let energyKineticR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
						let energyPotentialR = (xR - xmin) * externalForceX + (yR - ymin) * externalForceY;
						let energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
						let pressureR = (this.gamma - 1) * densityR * energyThermalR;
						let speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
						let hTotalR = energyTotalR + pressureR / densityR;
						let roeWeightR = Math.sqrt(densityR);

						let denom = roeWeightL + roeWeightR;
						let velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
						let velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
						let hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;
						
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
			let mindum = undefined;
			for (let j = 1; j < this.nx; ++j) {
				for (let i = 1; i < this.nx; ++i) {
					if (this.solid[i + this.nx * j]) continue;
					for (let side = 0; side < 2; ++side) {
						let maxLambda = Math.max(0, 
							this.interfaceEigenvalues[0+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[1+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[2+4*(side+2*(i+(this.nx+1)*j))],
							this.interfaceEigenvalues[3+4*(side+2*(i+(this.nx+1)*j))]);
						let minLambda = Math.min(0, 
							this.interfaceEigenvalues[0+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
							this.interfaceEigenvalues[1+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
							this.interfaceEigenvalues[2+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
							this.interfaceEigenvalues[3+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))]);
						let dx = this.xi[side + 2 * (i+dirs[side][0] + (this.nx+1) * (j+dirs[side][1]))] 
							- this.xi[side + 2 * (i + (this.nx+1) * j)];
						let dum = dx / (maxLambda - minLambda);
						if (mindum === undefined || dum < mindum) mindum = dum;
					}
				}
			}
			return this.cfl * mindum;
		}
	}
};

let HydroCell = function(args) {
	this.state = args.state;
	this.edges = args.edges;
	this.q = new Float64Array(4);	//state: rho, rho u, rho v, rho e
	this.pressure = 0;

	this.x = [0,0];
	for (let i = 0; i < this.edges.length; ++i) {
		let edge = this.state.edges[this.edges[i]];
		for (let j = 0; j < 2; ++j) {
			let vertex = this.state.vertexes[edge.vertexes[j]];
			this.x[0] += vertex[0];
			this.x[1] += vertex[1];
		}
	}
	this.x[0] /= 2 * this.edges.length;
	this.x[1] /= 2 * this.edges.length;
};

let HydroEdge = function(args) {
	this.state = args.state;
	this.vertexes = args.vertexes;
	
	let va = this.state.vertexes[this.vertexes[0]];
	let vb = this.state.vertexes[this.vertexes[1]];

	this.x = [
		.5 * (va[0] + vb[0]),
		.5 * (va[1] + vb[1])
	];

	let dx = vb[0] - va[0];
	let dy = vb[1] - va[1];
	//normalize
	let ds = Math.sqrt(dx * dx + dy * dy);
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
	for (let j = 0; j < 4; ++j) {
		this.eigenvalues[i] = 1;
		for (let i = 0; i < 4; ++i) {
			let delta_ij = i == j ? 1 : 0;
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
			let cellAX = this.state.cells[this.cells[0]].x;
			let cellBX = this.state.cells[this.cells[1]].x;
			let dax = [cellAX[0] - this.x[0], cellAX[1] - this.x[1]];
			let dbx = [cellBX[0] - this.x[0], cellBX[1] - this.x[1]];
			let aforward = dax[0] * this.normal[0] + dax[1] * this.normal[1] > 0;
			let bforward = dbx[0] * this.normal[0] + dbx[1] * this.normal[1] > 0;
			if (aforward && bforward) throw "both cells in front of normal";
			if (!aforward && !bforward) throw "both cells behind normal";
			if (!aforward && bforward) {
				let tmp = this.cells[0];
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

let HydroState = makeClass({ 
	init : function(args) {
		this.nx = args.size;
		this.cfl =.5;
		this.gamma = args.gamma;

		//define vertices
		
		this.vertexes = [];

		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				let x = i / this.nx * (xmax - xmin) + xmin;
				let y = j / this.nx * (ymax - ymin) + ymin;
				this.vertexes.push([x, y]);
			}
		}
		
		//define edges by vertices 
		
		this.edges = [];
		
		//x_{i-1/2},{j-1/2},dim: interface positions
		//0 <= i,j < this.nx+1
		//0 <= dim < 2
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					if ((side == 0 && i > 0) || (side == 1 && j > 0)) {
						let vertexes = [];
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
		let thiz = this;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				let edges = [];
				let cellIndex = this.cells.length;
				let addEdgeWithVertexes = function(a,b) {
					for (let edgeIndex = 0; edgeIndex < thiz.edges.length; ++edgeIndex) {
						let edge = thiz.edges[edgeIndex];
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
		for (let i = 0; i < this.edges.length; ++i) {
			this.edges[i].initAux();
		}
	
		this.resetSod();
		
		//qiTilde_{i-1/2},{j-1/2},side,state	
		this.interfaceDeltaQTilde = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4);
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					for (let state = 0; state < 4; ++state) {
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
		for (let i = 0; i < this.cells.length; ++i) {
			let cell = this.cells[i];
			let x = cell.x[0];
			let y = cell.x[1];
			let rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
			let u = 0;
			let v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			let energyKinetic = .5 * (u * u + v * v);
			let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			let energyThermal = 1;
			let energyTotal = energyKinetic + energyThermal + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u; 
			cell.q[2] = rho * v; 
			cell.q[3] = rho * energyTotal; 
			if (cell.q[0] !== cell.q[0]) throw 'nan';
		}
	},
	resetSodCylinder : function() {
		//TODO restructure the mesh to avoid the cylinder
		for (let i = 0; i < this.cells.length; ++i) {
			let cell = this.cells[i];
			let x = cell.x[0];
			let y = cell.x[1];
			let rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
			let u = 0;
			let v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			let energyKinetic = .5 * (u * u + v * v);
			let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			let energyThermal = 1;
			let energyTotal = energyKinetic + energyThermal + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u; 
			cell.q[2] = rho * v; 
			cell.q[3] = rho * energyTotal; 
		}
	},
	resetWave : function() {
		let xmid = .5 * (xmin + xmax);
		let ymid = .5 * (ymin + ymax);
		let dg = .2 * (xmax - xmin);
		for (let i = 0; i < this.cells.length; ++i) {
			let cell = this.cells[i];
			let x = cell.x[0];
			let y = cell.y[0];
			let dx = x - xmid;
			let dy = y - ymid;
			let rho = 3 * Math.exp(-(dx * dx + dy * dy) / (dg * dg)) + .1;
			let u = 0; 
			let v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			let energyKinetic = .5 * (u * u + v * v);
			let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			let energyThermal = 1;
			let energyTotal = energyKinetic + energyThermal + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u;
			cell.q[2] = rho * v;
			cell.q[3] = rho * energyTotal;
		}
	},
	//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
	resetKelvinHemholtz : function() {
		let xmid = .5 * (xmin + xmax);
		for (let i = 0; i < this.cells.length; ++i) {
			let cell = this.cells[i];
			let x = cell.x[0];
			let y = cell.x[1];
			let yInTheMiddle = y > (.75 * ymin + .25 * ymax) && y < (.25 * ymin + .75 * ymax);
			let rho = yInTheMiddle ? 2 : 1;
			let u = yInTheMiddle ? .5 : -.5;
			let v = 0;
			if (useNoise) {
				u += (Math.random() - .5) * 2 * noiseAmplitude;
				v += (Math.random() - .5) * 2 * noiseAmplitude;
			}
			//P = (gamma - 1) rho (eTotal - eKinetic - ePotential)
			//eTotal = P / ((gamma - 1) rho) + eKinetic + ePotential
			//eTotal = eThermal + eKinetic + ePotential
			//eThermal = P / ((gamma - 1) rho)
			let pressure = 2.5;
			let energyKinetic = .5 * (u * u + v * v);
			let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			let energyTotal = pressure / ((this.gamma - 1) * rho) + energyKinetic + energyPotential;
			cell.q[0] = rho;
			cell.q[1] = rho * u; 
			cell.q[2] = rho * v;
			cell.q[3] = rho * energyTotal; 
		}
		//TODO make it periodic on the left/right borders and reflecting on the top/bottom borders
	},
	//http://www.astro.virginia.edu/VITA/ATHENA/rt.html
	resetRayleighTaylor : function() {
		let ymid = .5 * (ymin + ymax);
		for (let i = 0; i < this.cells.length; ++i) {
			let cell = this.cells[i];
			let x = cell.x[0 + 2 * e];
			let y = cell.x[1 + 2 * e];
			let yGreaterThanMid = y > ymid;
			let rho = yGreaterThanMid ? 2 : 1;
			let u = 0;
			let v = 0;
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
			let energyPotential = (x - xmin) * externalForceX + (y - ymin) * externalForceY;
			let pressure = 2.5 - rho * energyPotential; 
			//let pressure = (this.gamma - 1) * rho * (2.5 - energyPotential);	//this fits with our non-external-force static ...
			/*
			pressure gradient zero ...
			dPressure / dy + rho * externalForce[side] = 0
			...integrate wrt y ...
			pressure + rho * externalForce.y * y = constant
			pressure = constant - rho * externalForce.y = constant - rho * energyPotential
			*/
			let energyKinetic = .5 * (u * u + v * v);
			/*
			energyTotal = energyKinetic + energyPotential + energyThermal
			energyTotal = .5 * (u*u + v*v) + energyPotential + pressure / (rho * (gamma - 1))
			*/
			let energyTotal = pressure / (rho * (this.gamma - 1)) + energyKinetic + energyPotential;
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
		let dt;
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


let Hydro = makeClass({
	init : function() {
		let size = Number($.url().param('size'));
		if (size === undefined || size !== size) size = 50;
		let gamma = Number($.url().param('gamma'));
		if (gamma === undefined || gamma !== gamma) gamma = 7/5;
		this.state = new HydroState({
			size : size,
			gamma : gamma
		});

		let plainShader = new GL.ShaderProgram({
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

		let stateShader = new GL.ShaderProgram({
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
		let vertexData = [];
		for (let i = 0; i < this.state.vertexes.length; ++i) {
			for (let j = 0; j < 2; ++j) {
				vertexData.push(this.state.vertexes[i][j]);
			}
		}
		let vertexBuffer = new GL.ArrayBuffer({dim : 2, data : vertexData, keep : true});
		this.vertexObj = new GL.SceneObject({
			mode : gl.POINTS,
			attrs : {
				vertex : vertexBuffer
			},
			shader : plainShader,
			static : true
		});
	
		//edges
		let indexData = [];
		for (let i = 0; i < this.state.edges.length; ++i) {
			for (let j = 0; j < 2; ++j) {
				indexData.push(this.state.edges[i].vertexes[j]);
			}
		}
		let indexBuffer = new GL.ElementArrayBuffer({data : indexData})
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
		let cellVertexData = [];
		let cellStateData = [];
		for (let i = 0; i < this.state.cells.length; ++i) {
			let cell = this.state.cells[i];
			for (let j = 0; j < 2; ++j) {
				cellVertexData.push(cell.x[j]);
			}
			cellStateData.push(0);
		}
		let cellVertexBuffer = new GL.ArrayBuffer({dim : 2, data : cellVertexData, keep : true});
		let cellStateBuffer = new GL.ArrayBuffer({dim : 1, data : cellStateData, keep : true});
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
		let result = drawToScreenMethods[this.state.drawToScreenMethod].call(this);
		let dataMin = result[0];
		let dataMax = result[1];
		
		if (this.updateLastDataRange) {
			this.lastDataMin = dataMin;
			this.lastDataMax = dataMax;
			let thisTick = Math.floor(Date.now() / 1000);
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

let hydro;
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
	let select = $('#' + id);
	for (let k in map) {
		let option = $('<option>', {text : k});
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

let colorSchemes = {};

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

	let panelContent = $('#content');
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
		let select = $('#gridsize');
		$.each([20, 50, 100, 200, 500, 1000], function(i,gridsize){
			let option = $('<option>', {text : gridsize});
			if (hydro.state.nx == gridsize) option.attr('selected', 'true');
			option.appendTo(select);
		});
		select.change(function(){
			let params = $.url().param();
			params.size = select.val();
			let url = location.href.match('[^?]*');
			let sep = '?';
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
			let v = Number($(this).val());
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
		let v = Number($(this).val());
		if (v != v) return;
		hydro.state.cfl = v;
	});
	$('#timeStepFixed').change(function() {
		if (!$(this).is(':checked')) return;
		useCFL = false;
	});
	$('#timeStepValue').val(fixedDT);
	$('#timeStepValue').change(function() {
		let v = Number($(this).val());
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

	let isobarSize = 16;
	let isobarData = new Uint8Array(isobarSize);
	for (let i = 1; i < isobarSize; i += 2) {
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

	for (let colorSchemeName in colorSchemes) {
		(function(){
			let v = colorSchemes[colorSchemeName];
			$('#colorScheme').append($('<option>', {
				value : colorSchemeName,
				text : colorSchemeName
			}));
		})();
	}
	$('#colorScheme').change(function() {
		let colorSchemeName = $(this).val();
		gl.bindTexture(gl.TEXTURE_2D, colorSchemes[colorSchemeName].obj);
	});

	let zoomFactor = .0003;	// upon mousewheel
	let dragging = false;
	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
			dragging = false;
		},
		move : function(dx,dy) {
			dragging = true;
			let aspectRatio = canvas.width / canvas.height;
			GL.view.pos[0] -= dx / canvas.width * 2 * (aspectRatio * GL.view.fovY);
			GL.view.pos[1] += dy / canvas.height * 2 * GL.view.fovY;
			GL.updateProjection();
		},
		zoom : function(zoomChange) {
			dragging = true;
			let scale = Math.exp(-zoomFactor * zoomChange);
			GL.view.fovY *= scale 
			GL.updateProjection();
		}
	});
	
	//start it off
	$(window).resize(onresize);
	onresize();
	update();
});
