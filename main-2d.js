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
import {DOM, getIDs, removeFromParent, show, hide, hidden} from '/js/util.js';
import {GLUtil} from '/js/gl-util.js';
import {makeGradient} from '/js/gl-util-Gradient.js';
import {Mouse3D} from '/js/mouse3d.js';
import {makeExplicitMethods} from './explicit_methods.js';
import {fluxMethods} from './flux_methods.js';

const ids = getIDs();
window.ids = ids;

const urlparams = new URLSearchParams(location.search);

let gl;
let glutil;
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

//global for storing stuff for key lookup
const _G = {};

_G.externalForceX = 0;
_G.externalForceY = 0;

_G.boundaryTopConstantValue = 0;
_G.boundaryLeftConstantValue = 0;
_G.boundaryRightConstantValue = 0;
_G.boundaryBottomConstantValue = 0;

let waveVtxBuf, waveStateBuf;

//interface directions
let dirs = [[1,0], [0,1]];

let drawToScreenMethods = {
	Density : function() {
		let nx = this.state.nx;
		let q = this.state.q;
		let dataMin = q[0];
		let dataMax = dataMin + 1e-9;
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		let e = 0;
		for (let j = 0; j < nx; ++j) {
			for (let i = 0; i < nx; ++i) {
				let s = q[0 + 4 * e];
				if (s < dataMin) dataMin = s;
				if (s > dataMax) dataMax = s;
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		return [dataMin, dataMax];
	},
	Velocity : function() {
		let nx = this.state.nx;
		let q = this.state.q;
		let dataMin = Math.sqrt(q[1] * q[1] + q[2] * q[2]) / q[0];
		let dataMax = dataMin + 1e-9; 
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		let e = 0;
		for (let j = 0; j < nx; ++j) {
			for (let i = 0; i < nx; ++i) {
				let rho = q[0 + 4 * e];
				let mx = q[1 + 4 * e];
				let my = q[2 + 4 * e];
				let s = Math.sqrt(mx * mx + my * my) / rho;
				if (s < dataMin) dataMin = s;
				if (s > dataMax) dataMax = s;
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		return [dataMin, dataMax];
	},
	Energy : function() {
		let nx = this.state.nx;
		let q = this.state.q;
		let dataMin = q[3] / q[0]; 
		let dataMax = dataMin + 1e-9; 
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		let e = 0;
		for (let j = 0; j < nx; ++j) {
			for (let i = 0; i < nx; ++i) {
				let s = q[3 + 4 * e] / q[0 + 4 * e]; 
				if (s < dataMin) dataMin = s;
				if (s > dataMax) dataMax = s;
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		return [dataMin, dataMax];
	},
	Pressure : function() {
		let nx = this.state.nx;
		let gamma = this.state.gamma;
		let q = this.state.q;
		let stateX = this.state.x;
		let x = stateX[0];
		let y = stateX[1];
		let rho = q[0];
		let u = q[1] / rho;
		let v = q[2] / rho;
		let energyTotal = q[3] / rho;
		let energyKinetic = .5 * (u * u + v * v);
		let energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
		let energyThermal = energyTotal - energyKinetic - energyPotential;
		let pressure = (gamma - 1) * rho * energyThermal;
		let dataMin = pressure;
		let dataMax = dataMin + 1e-9;
		let lastDataRange = this.lastDataMax - this.lastDataMin;
		for (let e = 1; e < nx * nx; ++e) {
			let x = stateX[0 + 2 * e];
			let y = stateX[1 + 2 * e];
			let rho = q[0 + 4 * e];
			let u = q[1 + 4 * e] / rho;
			let v = q[2 + 4 * e] / rho;
			let energyTotal = q[3 + 4 * e] / rho;
			let energyKinetic = .5 * (u * u + v * v);
			let energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
			let energyThermal = energyTotal - energyKinetic - energyPotential;
			let s = (gamma - 1) * rho * energyThermal;			
			if (s < dataMin) dataMin = s;
			if (s > dataMax) dataMax = s;
			this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
		}
		return [dataMin, dataMax];
	},
	Curl : function() {
		const nx = this.state.nx;
		const q = this.state.q;
		let dataMin = undefined; 
		let dataMax = undefined; 
		const lastDataRange = this.lastDataMax - this.lastDataMin;
		let e = 0;
		const dx = (xmax - xmin) / nx;
		const dy = (ymax - ymin) / nx;
		for (let j = 0; j < nx; ++j) {
			for (let i = 0; i < nx; ++i) {
				const rho = q[0 + 4 * e];
				const mx = q[1 + 4 * e];
				const my = q[2 + 4 * e];
				const mxPrev = this.state.q[1 + 4 * ( ((i+nx-1)%nx) + nx * j )];
				const myPrev = this.state.q[1 + 4 * ( i + nx * ((j+nx-1)%nx) )];
				const dmx = mx - mxPrev;
				const dmy = my - myPrev;
				const s = (dmx / dy - dmy / dx) / rho;
				if (dataMin === undefined || s < dataMin) dataMin = s;
				if (dataMax === undefined || s > dataMax) dataMax = s;
				this.vertexStates[e] = (s - this.lastDataMin) / lastDataRange;
				++e;
			}
		}
		return [dataMin, dataMax];
	}
};

const copyState = function(srcQ, destQ) {
	for (let i = 0; i < srcQ.length; ++i) {
		destQ[i] = srcQ[i];
	}
	return destQ;
};

const addMulState = function(to, from, scalar) {
	for (let i = 0; i < to.length; ++i) {
		to[i] += scalar * from[i];
	}
};

const explicitMethods = makeExplicitMethods({
	copyState : copyState,
	addMulState : addMulState,
});

let boundaryMethods = {
	//none : {top : function() {}, left : function() {}, right : function() {}, bottom : function() {}},
	periodic : {
		bottom : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					//top
					q[state + 4 * (i + nx * 0)] = q[state + 4 * (i + nx * (nx-4))];
					q[state + 4 * (i + nx * 1)] = q[state + 4 * (i + nx * (nx-3))];
				}
			}
		},
		left : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (0 + nx * i)] = q[state + 4 * (nx-4 + nx * i)];
					q[state + 4 * (1 + nx * i)] = q[state + 4 * (nx-3 + nx * i)];
				}
			}
		},
		right : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (nx-2 + nx * i)] = q[state + 4 * (2 + nx * i)];
					q[state + 4 * (nx-1 + nx * i)] = q[state + 4 * (3 + nx * i)];
				}
			}
		},
		top : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (i + nx * (nx-2))] = q[state + 4 * (i + nx * 2)];
					q[state + 4 * (i + nx * (nx-1))] = q[state + 4 * (i + nx * 3)];
				}
			}
		}
	},
	mirror : {
		bottom : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				q[0 + 4 * (i + nx * (0))] = q[0 + 4 * (i + nx * (3))];
				q[0 + 4 * (i + nx * (1))] = q[0 + 4 * (i + nx * (2))];
				q[1 + 4 * (i + nx * (0))] = -q[1 + 4 * (i + nx * (3))];
				q[1 + 4 * (i + nx * (1))] = -q[1 + 4 * (i + nx * (2))];
				q[2 + 4 * (i + nx * (0))] = -q[2 + 4 * (i + nx * (3))];
				q[2 + 4 * (i + nx * (1))] = -q[2 + 4 * (i + nx * (2))];
				q[3 + 4 * (i + nx * (0))] = q[3 + 4 * (i + nx * (3))];
				q[3 + 4 * (i + nx * (1))] = q[3 + 4 * (i + nx * (2))];
			}
		},
		left : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				q[0 + 4 * (0 + nx * i)] = q[0 + 4 * (3 + nx * i)];
				q[0 + 4 * (1 + nx * i)] = q[0 + 4 * (2 + nx * i)];
				q[1 + 4 * (0 + nx * i)] = -q[1 + 4 * (3 + nx * i)];
				q[1 + 4 * (1 + nx * i)] = -q[1 + 4 * (2 + nx * i)];
				q[2 + 4 * (0 + nx * i)] = -q[2 + 4 * (3 + nx * i)];
				q[2 + 4 * (1 + nx * i)] = -q[2 + 4 * (2 + nx * i)];
				q[3 + 4 * (0 + nx * i)] = q[3 + 4 * (3 + nx * i)];
				q[3 + 4 * (1 + nx * i)] = q[3 + 4 * (2 + nx * i)];
			}
		},
		right : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				q[0 + 4 * (nx-2 + nx * i)] = q[0 + 4 * (nx-3 + nx * i)];
				q[0 + 4 * (nx-1 + nx * i)] = q[0 + 4 * (nx-4 + nx * i)];
				q[1 + 4 * (nx-2 + nx * i)] = -q[1 + 4 * (nx-3 + nx * i)];
				q[1 + 4 * (nx-1 + nx * i)] = -q[1 + 4 * (nx-4 + nx * i)];
				q[2 + 4 * (nx-2 + nx * i)] = -q[2 + 4 * (nx-3 + nx * i)];
				q[2 + 4 * (nx-1 + nx * i)] = -q[2 + 4 * (nx-4 + nx * i)];
				q[3 + 4 * (nx-2 + nx * i)] = q[3 + 4 * (nx-3 + nx * i)];
				q[3 + 4 * (nx-1 + nx * i)] = q[3 + 4 * (nx-4 + nx * i)];
			}
		},
		top : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				q[0 + 4 * (i + nx * (nx-2))] = q[0 + 4 * (i + nx * (nx-3))];
				q[0 + 4 * (i + nx * (nx-1))] = q[0 + 4 * (i + nx * (nx-4))];
				q[1 + 4 * (i + nx * (nx-2))] = -q[1 + 4 * (i + nx * (nx-3))];
				q[1 + 4 * (i + nx * (nx-1))] = -q[1 + 4 * (i + nx * (nx-4))];
				q[2 + 4 * (i + nx * (nx-2))] = -q[2 + 4 * (i + nx * (nx-3))];
				q[2 + 4 * (i + nx * (nx-1))] = -q[2 + 4 * (i + nx * (nx-4))];
				q[3 + 4 * (i + nx * (nx-2))] = q[3 + 4 * (i + nx * (nx-3))];
				q[3 + 4 * (i + nx * (nx-1))] = q[3 + 4 * (i + nx * (nx-4))];
			}
		}
	},
	/*
	constant : {
		bottom : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state <= 2; ++state) {
					q[state + 4 * (i + nx * (0))] = q[state + 4 * (i + nx * (1))] = q[state + 4 * (i + nx * (2))];
				}
				for (let offset = 0; offset <= 1; ++offset) {
					q[3 + 4 * (i + nx * offset)] = _G.boundaryTopConstantValue;
				}
			}
		},
		left : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state <= 2; ++state) {
					q[state + 4 * (0 + nx * i)] = q[state + 4 * (1 + nx * i)] = q[state + 4 * (2 + nx * i)];
				}
				for (let offset = 0; offset <= 1; ++offset) {
					q[3 + 4 * (offset + nx * i)] = _G.boundaryLeftConstantValue;
				}
			}
		},
		right : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state <= 2; ++state) {
					q[state + 4 * (nx-1 + nx * i)] = q[state + 4 * (nx-2 + nx * i)] = q[state + 4 * (nx-3 + nx * i)];
				}
				for (let offset = 0; offset <= 1; ++offset) {
					q[3 + 4 * (nx-1-offset + nx * i)] = _G.boundaryRightConstantValue;
				}
			}
		},
		top : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state <= 2; ++state) {
					q[state + 4 * (i + nx * (nx-1))] = q[state + 4 * (i + nx * (nx-2))] = q[state + 4 * (i + nx * (nx-3))];
				}
				for (let offset = 0; offset <= 1; ++offset) {
					q[3 + 4 * (i + nx * (nx-1-offset))] = _G.boundaryBottomConstantValue;
				}
			}
		}
	},
	*/
	freeflow : {
		bottom : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (i + nx * (0))] = q[state + 4 * (i + nx * (1))] = q[state + 4 * (i + nx * (2))];
				}
			}
		},
		left : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (0 + nx * i)] = q[state + 4 * (1 + nx * i)] = q[state + 4 * (2 + nx * i)];
				}
			}
		},
		right : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (nx-1 + nx * i)] = q[state + 4 * (nx-2 + nx * i)] = q[state + 4 * (nx-3 + nx * i)];
				}
			}
		},
		top : function(nx,q) {
			for (let i = 0; i < nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					q[state + 4 * (i + nx * (nx-1))] = q[state + 4 * (i + nx * (nx-2))] = q[state + 4 * (i + nx * (nx-3))];
				}
			}
		}
	}
};

class EulerEquationBurgersSolver {
	initStep() {}
	calcCFLTimestep() {
		let mindum = undefined;
		let qIndex = 0;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				if (!this.solid[i + this.nx * j]) { 
					const x = this.x[0 + 2 * (i + this.nx * j)];
					const y = this.x[1 + 2 * (i + this.nx * j)];
					const rho = this.q[0 + qIndex];
					const u = this.q[1 + qIndex] / rho;
					const v = this.q[2 + qIndex] / rho; 
					const energyTotal = this.q[3 + qIndex] / rho; 
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = energyTotal - energyKinetic - energyPotential;
					const speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyThermal);
					const dx = this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] - this.xi[0 + 2 * (i + (this.nx+1) * j)];
					const dy = this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] - this.xi[1 + 2 * (i + (this.nx+1) * j)];
					let dum = dx / (speedOfSound + Math.abs(u));
					if (mindum === undefined || dum < mindum) mindum = dum;
					dum = dy / (speedOfSound + Math.abs(v));
					if (mindum === undefined || dum < mindum) mindum = dum;
				}
				qIndex += 4;
			}
		}
		return this.cfl * mindum;
	}
}

class EulerEquationBurgersExplicit extends EulerEquationBurgersSolver {
	step(dt) {
		//get velocity at interfaces from state
		for (let j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
			for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
				let uiIndex = 2 * (i + (this.nx+1) * j);
				for (let side = 0; side < 2; ++side) {
					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					const indexR = i + this.nx * j;
					const qIndexL = 4 * indexL;
					const qIndexR = 4 * indexR;
					const solidL = this.solid[indexL];
					const solidR = this.solid[indexR];
					let uL = this.q[1+side + qIndexL] / this.q[0 + qIndexL];
					let uR = this.q[1+side + qIndexR] / this.q[0 + qIndexR];
					if (solidL && !solidR) {
						uL = -uR;
					} else if (solidR && !solidL) {
						uR = -uL;
					}
					this.ui[side + uiIndex] = .5 * (uL + uR);
				}
			}
		}
		//boundary zero
		for (let j = 0; j < this.nghost; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					//left boundary, left and top sides, zero vector
					this.ui[side + 2 * (j + (this.nx+1) * i)] = 0;
					//right boundary, left and top sides, zero vector
					this.ui[side + 2 * (this.nx-j + (this.nx+1) * i)] = 0;
					//top boundary, left and top sides, zero vector
					this.ui[side + 2 * (i + (this.nx+1) * j)] = 0;
					//bottom boundary, left and top sides, zero vector
					this.ui[side + 2 * (i + (this.nx+1) * (this.nx-j))] = 0;
				}
			}
		}

		//compute flux and advect for each state vector
		for (let state = 0; state < 4; ++state) {	//state
			//r_{i-1/2},{j-1/2} flux limiter
			for (let j = this.nghost; j < this.nx+this.nghost-3; ++j) {
				for (let i = this.nghost; i < this.nx+this.nghost-3; ++i) {
					for (let side = 0; side < 2; ++side) {
						const indexL2 = i - 2*dirs[side][0] + this.nx * (j - 2*dirs[side][1]);
						const indexL1 = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
						const indexR1 = i + this.nx * j;
						const indexR2 = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
						
						const solidL2 = this.solid[indexL2];
						const solidL1 = this.solid[indexL1];
						const solidR1 = this.solid[indexR1];
						const solidR2 = this.solid[indexR2];
						
						let qL2 = this.q[state + 4 * indexL2];
						let qL1 = this.q[state + 4 * indexL1];
						let qR1 = this.q[state + 4 * indexR1];
						let qR2 = this.q[state + 4 * indexR2];
						
						/*
						slope limiters and arbitrary boundaries...
						we have to look two over in any direction
						if there's a wall two over, what do we do?
						*/
						if (solidL2) {
							qL2 = qL1;
							if (state == 1+side) qL2 = -qL2;
						}
						if (solidR2) {
							qR2 = qR1;
							if (state == 1+side) qR2 = -qR2;
						}
						if (solidL1) {
							qL1 = -qR1;
							if (state == 1+side) qL1 = -qL1;
						}
						if (solidR1) {
							qR1 = -qL1;
							if (state == 1+side) qR1 = -qR1;
						}

						//dq = q_i,j - q_{{i,j}-dirs[side]}
						let dq = qR1 - qL1;
						if (Math.abs(dq) > 0) {
							if (this.ui[side + 2 * (i + (this.nx+1) * j)] >= 0) {
								this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = (qL1 - qL2) / dq;
							} else {
								this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = (qR2 - qR1) / dq;
							}
						} else {
							this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
						}
					}
				}
			}
		
			//now for ghost boundaries
			for (let j = 0; j < this.nghost; ++j) {
				for (let i = 0; i <= this.nx; ++i) {
					for (let side = 0; side < 2; ++side) {
						this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
						this.r[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
						this.r[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
						this.r[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
					}
				}
			}


			let dxi = [];
			//construct flux:
			for (let j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
				for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
					const dx = this.x[0 + 2 * (i + this.nx * j)] 
						- this.x[0 + 2 * (i-dirs[0][0] + this.nx * (j-dirs[0][1]))];
					const dy = this.x[1 + 2 * (i + this.nx * j)] 
						- this.x[1 + 2 * (i-dirs[1][0] + this.nx * (j-dirs[1][1]))];
					dxi[0] = dx;
					dxi[1] = dy;
					const volume = dx * dy;
					//flux calculation 
					for (let side = 0; side < 2; ++side) {
						
						const rIndex = state + 4 * (side + 2 * (i + (this.nx+1) * j));
						const uiIndex = side + 2 * (i + (this.nx+1) * j);
						const fluxIndex = state + 4 * uiIndex;
						const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
						const indexR = i + this.nx * j;
						const qIndexL = state + 4 * indexL;
						const qIndexR = state + 4 * indexR;
						
						//apply limiter
						const phi = fluxMethods[this.fluxMethod](this.r[rIndex]);
						const ui = this.ui[uiIndex];
						
						const solidL = this.solid[indexL];
						const solidR = this.solid[indexR];
						
						
						let qL = this.q[qIndexL];
						let qR = this.q[qIndexR];
						if (solidL && !solidR) {
							qL = qR;
							if (state == side+1) qL = -qL;
						}
						if (solidR && !solidL) {
							qR = qL;
							if (state == side+1) qR = -qR;
						}
						if (ui >= 0) {
							this.flux[fluxIndex] = ui * qL; 
						} else {
							this.flux[fluxIndex] = ui * qR; 
						}
						let delta = phi * (qR - qL);
						this.flux[fluxIndex] += delta * .5 * .5
							* Math.abs(ui)
							* (1 - Math.abs(ui * dt / dxi[side]));
					}
				}
			}
			
			//now for ghost boundaries
			for (let j = 0; j < this.nghost-1; ++j) {
				for (let i = 0; i <= this.nx; ++i) {
					for (let side = 0; side < 2; ++side) {
						this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
						this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
						this.flux[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
						this.flux[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
					}
				}
			}

			//update cells
			dxi = [];
			for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
				for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
					if (this.solid[i + this.nx * j]) continue;

					let xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
					let xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
					const dx = this.xi[xiIndexR] - this.xi[xiIndexL];
					
					xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
					xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
					const dy = this.xi[xiIndexR] - this.xi[xiIndexL];
				
					const volume = dx * dy;
					dxi[0] = dx;
					dxi[1] = dy;
					
					for (let side = 0; side < 2; ++side) {
						const ifluxR = state + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])));
						const ifluxL = state + 4 * (side + 2 * (i + (this.nx+1) * j));
						const df = this.flux[ifluxR] - this.flux[ifluxL];
						this.q[state + 4 * (i + this.nx * j)] -= dt * df / dxi[side];//* volume / (dxi[side] * dxi[side]);
					}
				}
			}
		}
	
		//boundary again
		this.boundary();

		//compute pressure
		let qIndex = 0;
		let index = 0;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				if (!this.solid[index]) {
					const x = this.x[0 + 2 * (i + this.nx * j)];
					const y = this.x[1 + 2 * (i + this.nx * j)];
					const rho = this.q[0 + qIndex];
					const u = this.q[1 + qIndex] / rho; 
					const v = this.q[2 + qIndex] / rho; 
					const energyTotal = this.q[3 + qIndex] / rho; 
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = energyTotal - energyKinetic - energyPotential;
					this.pressure[index] = (this.gamma - 1) * rho * energyThermal;
				}
				++index;
				qIndex += 4;
			}
		}

		//apply external force
		for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				const index = i + this.nx * j;
				const qIndex = 4 * index;
				const rho = this.q[0 + qIndex];
				this.q[3 + qIndex] -= dt * (_G.externalForceX * this.q[1 + qIndex] + _G.externalForceY * this.q[2 + qIndex]);
				this.q[1 + qIndex] -= dt * rho * _G.externalForceX;
				this.q[2 + qIndex] -= dt * rho * _G.externalForceY;
			}
		}

		//apply momentum diffusion = pressure
		for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				let index = i + this.nx * j;
				let qIndex = 4 * index;
				if (this.solid[index]) continue;
				for (let side = 0; side < 2; ++side) {
					const indexR = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					
					const pressureC = this.pressure[index];
					let pressureL = this.pressure[indexL];
					if (this.solid[indexL]) pressureL = pressureC;
					let pressureR = this.pressure[indexR];
					if (this.solid[indexR]) pressureR = pressureC;
					
					const dPressure = pressureR - pressureL;
					const dx = this.x[side + 2 * indexR] - this.x[side + 2 * indexL];
					const rho = this.q[0 + qIndex];
					this.q[1+side + qIndex] -= dt * dPressure / dx;
				}
			}
		}

		//apply work diffusion = momentum
		let dxi = [];
		for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				let index = i + this.nx * j;
				let qIndex = 4 * index;
				if (this.solid[index]) continue;
				for (let side = 0; side < 2; ++side) {
					const indexR = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
				
					const uC = this.q[1+side + 4 * index] / this.q[0 + 4 * index];

					//this is pulling the coordinate associated with the interface's direction
					//a more robust method would be to take both velocity components and dot them with the interface tangent
					let uR = this.q[1+side + 4 * indexR] / this.q[0 + 4 * indexR];
					if (this.solid[indexR]) {
						uR = -uC;
					}
					
					let uL = this.q[1+side + 4 * indexL] / this.q[0 + 4 * indexL];
					if (this.solid[indexL]) {
						uL = -uC;
					}
					
					const pressureC = this.pressure[index];
					let pressureL = this.pressure[indexL];
					if (this.solid[indexL]) pressureL = pressureC;
					let pressureR = this.pressure[indexR];
					if (this.solid[indexR]) pressureR = pressureC;
					
					let dx = this.x[side + 2 * indexR] - this.x[side + 2 * indexL];
					this.q[3 + qIndex] -= dt * (pressureR * uR - pressureL * uL) / dx;
				}
			}
		}
	}
}

class EulerEquationBurgersBackwardEulerGaussSeidel extends EulerEquationBurgersSolver {
	step(dt) {
		this.oldQ.set(this.q);
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			//get velocity at interfaces from state
			for (let j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
				for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
					const uiIndex = 2 * (i + (this.nx+1) * j);
					for (let side = 0; side < 2; ++side) {
						const qIndexR = 4 * (i + this.nx * j);
						const qIndexL = 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
						this.ui[side + uiIndex] = .5 * (
							this.q[1+side + qIndexR] / this.q[0 + qIndexR] + 
							this.q[1+side + qIndexL] / this.q[0 + qIndexL]);
					}
				}
			}
			//boundary zero
			for (let j = 0; j < this.nghost; ++j) {
				for (let i = 0; i <= this.nx; ++i) {
					for (let side = 0; side < 2; ++side) {
						//left boundary, left and top sides, zero vector
						this.ui[side + 2 * (j + (this.nx+1) * i)] = 0;
						//right boundary, left and top sides, zero vector
						this.ui[side + 2 * (this.nx-j + (this.nx+1) * i)] = 0;
						//top boundary, left and top sides, zero vector
						this.ui[side + 2 * (i + (this.nx+1) * j)] = 0;
						//bottom boundary, left and top sides, zero vector
						this.ui[side + 2 * (i + (this.nx+1) * (this.nx-j))] = 0;
					}
				}
			}

			//compute flux and advect for each state vector
			for (let state = 0; state < 4; ++state) {	//state
				//r_{i-1/2},{j-1/2} flux limiter
				for (let j = this.nghost; j < this.nx+this.nghost-3; ++j) {
					for (let i = this.nghost; i < this.nx+this.nghost-3; ++i) {
						for (let side = 0; side < 2; ++side) {
							//dq = q_i,j - q_{{i,j}-dirs[side]}
							let dq = this.q[state + 4 * (i + this.nx * j)]
								- this.q[state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))];
							if (Math.abs(dq) > 0) {
								if (this.ui[side + 2 * (i + (this.nx+1) * j)] >= 0) {
									this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										(this.q[state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]))]
											- this.q[state + 4 * (i - 2*dirs[side][0] + this.nx * (j - 2*dirs[side][1]))]) / dq;
								} else {
									this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 
										(this.q[state + 4 * (i + dirs[side][0] + this.nx * (j + dirs[side][1]))]
											- this.q[state + 4 * (i + this.nx * j)]) / dq;
									
								}
							} else {
								this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							}
						}
					}
				}
			
				//now for ghost boundaries
				for (let j = 0; j < this.nghost; ++j) {
					for (let i = 0; i <= this.nx; ++i) {
						for (let side = 0; side < 2; ++side) {
							this.r[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.r[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.r[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.r[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						}
					}
				}


				let dxi = [];
				//construct flux:
				for (let j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
					for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
						const dx = this.x[0 + 2 * (i + this.nx * j)] 
							- this.x[0 + 2 * (i-dirs[0][0] + this.nx * (j-dirs[0][1]))];
						const dy = this.x[1 + 2 * (i + this.nx * j)] 
							- this.x[1 + 2 * (i-dirs[1][0] + this.nx * (j-dirs[1][1]))];
						dxi[0] = dx;
						dxi[1] = dy;
						const volume = dx * dy;
						//flux calculation 
						for (let side = 0; side < 2; ++side) {
							
							const rIndex = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							//apply limiter
							const phi = fluxMethods[this.fluxMethod](this.r[rIndex]);
							const uiIndex = side + 2 * (i + (this.nx+1) * j);
							const fluxIndex = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							const qIndexL = state + 4 * (i - dirs[side][0] + this.nx * (j - dirs[side][1]));
							const qIndexR = state + 4 * (i + this.nx * j);
							if (this.ui[uiIndex] >= 0) {
								this.flux[fluxIndex] = this.ui[uiIndex] * this.q[qIndexL];
							} else {
								this.flux[fluxIndex] = this.ui[uiIndex] * this.q[qIndexR];
							}
							const delta = phi * (this.q[qIndexR] - this.q[qIndexL]);
							this.flux[fluxIndex] += delta * .5 * .5
								* Math.abs(this.ui[uiIndex])
								* (1 - Math.abs(this.ui[uiIndex] * dt / dxi[side]));//* volume / (dxi[side] * dxi[side])));
						}
					}
				}
				
				//now for ghost boundaries
				for (let j = 0; j < this.nghost-1; ++j) {
					for (let i = 0; i <= this.nx; ++i) {
						for (let side = 0; side < 2; ++side) {
							this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
							this.flux[state + 4 * (side + 2 * (i + (this.nx+1) * (this.nx-j)))] = 0;
							this.flux[state + 4 * (side + 2 * (j + (this.nx+1) * i))] = 0;
							this.flux[state + 4 * (side + 2 * (this.nx-j + (this.nx+1) * i))] = 0;
						}
					}
				}

				//update cells
				dxi = [];
				for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
					for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
						
						let xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
						let xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
						const dx = this.xi[xiIndexR] - this.xi[xiIndexL];
						
						xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
						xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
						const dy = this.xi[xiIndexR] - this.xi[xiIndexL];
					
						const volume = dx * dy;
						dxi[0] = dx;
						dxi[1] = dy;
						
						for (let side = 0; side < 2; ++side) {
							const ifluxR = state + 4 * (side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])));
							const ifluxL = state + 4 * (side + 2 * (i + (this.nx+1) * j));
							const df = this.flux[ifluxR] - this.flux[ifluxL];
							this.q[state + 4 * (i + this.nx * j)] = 
								this.oldQ[state + 4 * (i + this.nx * j)] 
								- dt * df / dxi[side];//* volume / (dxi[side] * dxi[side]);
						}
					}
				}
			}
		}
		
		//boundary again
		this.boundary();

		//apply external force
		this.oldQ.set(this.q);
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
				for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
					const qIndex = 4 * (i + this.nx * j);
					const rho = this.q[0 + qIndex];
					this.q[3 + qIndex] = this.oldQ[3 + qIndex] - dt * (_G.externalForceX * this.q[1 + qIndex] + _G.externalForceY * this.q[2 + qIndex]);
					this.q[1 + qIndex] = this.oldQ[1 + qIndex] - dt * rho * _G.externalForceX;
					this.q[2 + qIndex] = this.oldQ[2 + qIndex] - dt * rho * _G.externalForceY;
				}
			}
		}

		//apply momentum diffusion = pressure
		this.oldQ.set(this.q);
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			let qIndex = 0;
			let pIndex = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * (i + this.nx * j)];
					const y = this.x[1 + 2 * (i + this.nx * j)];
					const rho = this.q[0 + qIndex];
					const u = this.q[1 + qIndex] / rho; 
					const v = this.q[2 + qIndex] / rho; 
					const energyTotal = this.q[3 + qIndex] / rho; 
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = energyTotal - energyKinetic - energyPotential;
					this.pressure[pIndex] = (this.gamma - 1) * rho * energyThermal;
					++pIndex;
					qIndex += 4;
				}
			}			
		
			for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
				for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
					const qIndex = 4 * (i + this.nx * j);
					for (let side = 0; side < 2; ++side) {
						const plusIndex = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
						const minusIndex = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
						const dPressure = this.pressure[plusIndex] - this.pressure[minusIndex];
						const dx = this.x[side + 2 * plusIndex] - this.x[side + 2 * minusIndex];
						this.q[1+side + qIndex] = this.oldQ[1+side + qIndex] - dt * dPressure / dx;
					}
				}
			}
		}

		//apply work diffusion = momentum
		this.oldQ.set(this.q);
		let externalForce = [_G.externalForceX, _G.externalForceY];
		let dxi = [];
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			let qIndex = 0;
			let pIndex = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * (i + this.nx * j)];
					const y = this.x[1 + 2 * (i + this.nx * j)];
					const rho = this.q[0 + qIndex];
					const u = this.q[1 + qIndex] / rho; 
					const v = this.q[2 + qIndex] / rho; 
					const energyTotal = this.q[3 + qIndex] / rho; 
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = energyTotal - energyKinetic - energyPotential;
					this.pressure[pIndex] = (this.gamma - 1) * rho * energyThermal;
					++pIndex;
					qIndex += 4;
				}
			}			
			for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
				for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
					let qIndex = 4 * (i + this.nx * j);
					for (let side = 0; side < 2; ++side) {
						const plusIndex = i + dirs[side][0] + this.nx * (j + dirs[side][1]);
						const minusIndex = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
						const dx = this.x[side + 2 * plusIndex] - this.x[side + 2 * minusIndex];
						//this is pulling the coordinate associated with the interface's direction
						//a more robust method would be to take both velocity components and dot them with the interface tangent
						const uR = this.q[1+side + 4 * plusIndex] / this.q[0 + 4 * plusIndex];
						const uL = this.q[1+side + 4 * minusIndex] / this.q[0 + 4 * minusIndex];
						const pressureR = this.pressure[plusIndex];
						const pressureL = this.pressure[minusIndex];
						this.q[3 + qIndex] = this.oldQ[3 + qIndex] - dt * (pressureR * uR - pressureL * uL) / dx;
					}
				}
			}		
		}
	}
}

class GodunovSolver {
	calcCFLTimestep() {
		let mindum = undefined;
		for (let j = 1; j < this.nx; ++j) {
			for (let i = 1; i < this.nx; ++i) {
				if (this.solid[i + this.nx * j]) continue;
				for (let side = 0; side < 2; ++side) {
					const maxLambda = Math.max(0, 
						this.interfaceEigenvalues[0+4*(side+2*(i+(this.nx+1)*j))],
						this.interfaceEigenvalues[1+4*(side+2*(i+(this.nx+1)*j))],
						this.interfaceEigenvalues[2+4*(side+2*(i+(this.nx+1)*j))],
						this.interfaceEigenvalues[3+4*(side+2*(i+(this.nx+1)*j))]);
					const minLambda = Math.min(0, 
						this.interfaceEigenvalues[0+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
						this.interfaceEigenvalues[1+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
						this.interfaceEigenvalues[2+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))],
						this.interfaceEigenvalues[3+4*(side+2*(i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1])))]);
					const dx = this.xi[side + 2 * (i+dirs[side][0] + (this.nx+1) * (j+dirs[side][1]))] 
						- this.xi[side + 2 * (i + (this.nx+1) * j)];
					const dum = dx / (maxLambda - minLambda);
					if (mindum === undefined || dum < mindum) mindum = dum;
				}
			}
		}
		return this.cfl * mindum;
	}
	step(dt) {
		const deriv = GodunovSolver.prototype.calcDerivative;
		explicitMethods[this.explicitMethod].call(this, dt, deriv);
	}
	calcDerivative(dt, dq_dt) {
		const qLs = [];
		const qRs = [];
		for (let j = 1; j < this.nx; ++j) {
			for (let i = 1; i < this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					const indexR = i + this.nx * j;
					const interfaceIndex = side + 2 * (i + (this.nx+1) * j);
					const solidL = this.solid[indexL];
					const solidR = this.solid[indexR];
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
					const indexL2 = i - 2 * dirs[side][0] + this.nx * (j - 2 * dirs[side][1]);
					const indexL1 = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					const indexR1 = i + this.nx * j;
					const indexR2 = i + dirs[side][0] + this.nx * (j + dirs[side][1]); 
					
					const solidL2 = this.solid[indexL2];
					const solidL1 = this.solid[indexL1];
					const solidR1 = this.solid[indexR1];
					const solidR2 = this.solid[indexR2];
					
					const interfaceIndexL = side + 2 * (i - dirs[side][0] + (this.nx+1) * (j - dirs[side][1]));
					const interfaceIndex = side + 2 * (i + (this.nx+1) * j);
					const interfaceIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
					
					for (let state = 0; state < 4; ++state) {
						
						let interfaceDeltaQTildeL = this.interfaceDeltaQTilde[state + 4 * interfaceIndexL];
						const interfaceDeltaQTilde = this.interfaceDeltaQTilde[state + 4 * interfaceIndex];
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
	
		const fluxAvg = [];	//4
		const fluxTilde = [];	//4
		const dxi = [];
		//transform cell q's into cell qTilde's (eigenspace)
		for (let j = this.nghost-1; j < this.nx+this.nghost-2; ++j) {
			for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
				const dx = this.xi[0 + 2 * (i + (this.nx+1) * j)] 
					- this.xi[0 + 2 * (i-dirs[0][0] + (this.nx+1) * (j-dirs[0][1]))];
				const dy = this.xi[1 + 2 * (i + (this.nx+1) * j)] 
					- this.xi[1 + 2 * (i-dirs[1][0] + (this.nx+1) * (j-dirs[1][1]))];
				const volume = dx * dy;
				dxi[0] = dx;
				dxi[1] = dy;
				for (let side = 0; side < 2; ++side) {

					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					const indexR = i + this.nx * j;
					const solidL = this.solid[indexL];
					const solidR = this.solid[indexR];
					const interfaceIndex = side + 2 * (i + (this.nx+1) * j);
					
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
						const eigenvalue = this.interfaceEigenvalues[state + 4 * interfaceIndex];
						if (eigenvalue >= 0) {
							theta = 1;
						} else {
							theta = -1;
						}
					
						const phi = fluxMethods[this.fluxMethod](this.rTilde[state + 4 * interfaceIndex]);
						const epsilon = eigenvalue * dt / dxi[side];//* volume / (dxi[side] * dxi[side]); 
						const deltaFluxTilde = eigenvalue * this.interfaceDeltaQTilde[state + 4 * interfaceIndex];
						fluxTilde[state] = -.5 * deltaFluxTilde * (theta + .5 * phi * (epsilon - theta));
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
				
		for (let j = 0; j < this.nghost; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					//deriv while we're here
					//left
					dq_dt[state + 4 * (j + this.nx * i)] = 0;
					//right
					dq_dt[state + 4 * (this.nx-1-j + this.nx * i)] = 0;
					//bottom
					dq_dt[state + 4 * (i + this.nx * j)] = 0;
					//top
					dq_dt[state + 4 * (i + this.nx * (this.nx-1-j))] = 0;
				}
			}
		}

		//update cells
		for (let j = this.nghost; j < this.nx-this.nghost; ++j) {
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				for (let state = 0; state < 4; ++state) {
					dq_dt[state + 4 * (i + this.nx * j)] = 0;
				}
				
				if (this.solid[i + this.nx * j]) {
					continue;
				}

				let xiIndexR = 0 + 2 * (i + dirs[0][0] + (this.nx+1) * (j + dirs[0][1]));
				let xiIndexL = 0 + 2 * (i + (this.nx+1) * j);
				const dx = this.xi[xiIndexR] - this.xi[xiIndexL];
				
				xiIndexR = 1 + 2 * (i + dirs[1][0] + (this.nx+1) * (j + dirs[1][1]));
				xiIndexL = 1 + 2 * (i + (this.nx+1) * j);
				const dy = this.xi[xiIndexR] - this.xi[xiIndexL];
				
				const volume = dx * dy;
				dxi[0] = dx;
				dxi[1] = dy;
				
				for (let side = 0; side < 2; ++side) {
					const interfaceIndexR = side + 2 * (i + dirs[side][0] + (this.nx+1) * (j + dirs[side][1]));
					const interfaceIndexL = side + 2 * (i + (this.nx+1) * j);
					for (let state = 0; state < 4; ++state) {
						const interfaceFluxIndexR = state + 4 * interfaceIndexR;
						const interfaceFluxIndexL = state + 4 * interfaceIndexL;
						const df = this.flux[interfaceFluxIndexR] - this.flux[interfaceFluxIndexL];
						dq_dt[state + 4 * (i + this.nx * j)] -= df / dxi[side];//* volume / (dxi[side] * dxi[side]);
					}
				}
			}
		}
	}
}

class EulerEquationGodunovSolver extends GodunovSolver {
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
	buildEigenstate(offset, matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocityX, velocityY, hTotal, gamma, normalX, normalY) {

		if ((hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)) < 0) {
			console.log('sqrt error');
		}

		//calculate matrix & eigenvalues & vectors at interface from state at interface
		const speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * (velocityX * velocityX + velocityY * velocityY)));
		const tangentX = -normalY;
		const tangentY = normalX;
		const velocityN = velocityX * normalX + velocityY * normalY;
		const velocityT = velocityX * tangentX + velocityY * tangentY;
		const velocitySq = velocityX * velocityX + velocityY * velocityY;	
		
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
				const epsilon = 1e-5;
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
}

class EulerEquationGodunovExplicit extends EulerEquationGodunovSolver {
	initStep() {
		for (let j = 1; j < this.nx; ++j) {
			for (let i = 1; i < this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					const xL = side == 0 ? this.xi[0 + 2 * (i + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
					const xR = side == 0 ? this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
					const yL = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * j)] : this.x[1 + 2 * (i + this.nx * j)];
					const yR = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] : this.x[1 + 2 * (i + this.nx * j)];
			
					const normalX = dirs[side][0];
					const normalY = dirs[side][1];

					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					const qIndexL = 4 * indexL;
					const solidL = this.solid[indexL];
					const densityL = this.q[0 + qIndexL];
					const velocityXL = this.q[1 + qIndexL] / densityL;
					const velocityYL = this.q[2 + qIndexL] / densityL;
					const energyTotalL = this.q[3 + qIndexL] / densityL;
					
					const indexR = i + this.nx * j;
					const qIndexR = 4 * indexR;
					const solidR = this.solid[indexR];
					const densityR = this.q[0 + qIndexR];
					const velocityXR = this.q[1 + qIndexR] / densityR;
					const velocityYR = this.q[2 + qIndexR] / densityR;
					const energyTotalR = this.q[3 + qIndexR] / densityR;					
					
					if (solidL && !solidR) {	//right cell has a wall on the left
						densityL = densityR;
						//reflect along normal
						const velocityNR = normalX * velocityXR + normalY * velocityYR;
						velocityXL = velocityXR - 2 * normalX * velocityNR;
						velocityYL = velocityYR - 2 * normalY * velocityNR;
						energyTotalL = energyTotalR;
					} else if (!solidL && solidR) {	//left cell has a wall on the right
						densityR = densityL;
						//reflect along normal
						const velocityNL = normalX * velocityXL + normalY * velocityYL;
						velocityXR = velocityXL - 2 * normalX * velocityNL;
						velocityYR = velocityYL - 2 * normalY * velocityNL;
						energyTotalR = energyTotalL;
					}
				
					const energyKineticL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
					const energyPotentialL = (xL - xmin) * _G.externalForceX + (yL - ymin) * _G.externalForceY;
					const energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
					const pressureL = (this.gamma - 1) * densityL * energyThermalL;
					const speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
					const hTotalL = energyTotalL + pressureL / densityL;
					const roeWeightL = Math.sqrt(densityL);
				
					const energyKineticR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
					const energyPotentialR = (xR - xmin) * _G.externalForceX + (yR - ymin) * _G.externalForceY;
					const energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
					const pressureR = (this.gamma - 1) * densityR * energyThermalR;
					const speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
					const hTotalR = energyTotalR + pressureR / densityR;
					const roeWeightR = Math.sqrt(densityR);

					const denom = roeWeightL + roeWeightR;
					const velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
					const velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
					const hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;
					
					EulerEquationGodunovSolver.prototype.buildEigenstate(
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
	}
}

class EulerEquationRoeExplicit extends EulerEquationGodunovSolver {
	initStep() {
		for (let j = 1; j < this.nx; ++j) {
			for (let i = 1; i < this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					const xL = side == 0 ? this.xi[0 + 2 * (i + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
					const xR = side == 0 ? this.xi[0 + 2 * (i+1 + (this.nx+1) * j)] : this.x[0 + 2 * (i + this.nx * j)];
					const yL = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * j)] : this.x[1 + 2 * (i + this.nx * j)];
					const yR = side == 1 ? this.xi[1 + 2 * (i + (this.nx+1) * (j+1))] : this.x[1 + 2 * (i + this.nx * j)];
			
					const normalX = dirs[side][0];
					const normalY = dirs[side][1];

					const indexL = i - dirs[side][0] + this.nx * (j - dirs[side][1]);
					const qIndexL = 4 * indexL;
					const solidL = this.solid[indexL];
					let densityL = this.q[0 + qIndexL];
					let velocityXL = this.q[1 + qIndexL] / densityL;
					let velocityYL = this.q[2 + qIndexL] / densityL;
					let energyTotalL = this.q[3 + qIndexL] / densityL;
					
					const indexR = i + this.nx * j;
					const qIndexR = 4 * indexR;
					const solidR = this.solid[indexR];
					let densityR = this.q[0 + qIndexR];
					let velocityXR = this.q[1 + qIndexR] / densityR;
					let velocityYR = this.q[2 + qIndexR] / densityR;
					let energyTotalR = this.q[3 + qIndexR] / densityR;					
					
					if (solidL && !solidR) {	//right cell has a wall on the left
						densityL = densityR;
						//reflect along normal
						const velocityNR = normalX * velocityXR + normalY * velocityYR;
						velocityXL = velocityXR - 2 * normalX * velocityNR;
						velocityYL = velocityYR - 2 * normalY * velocityNR;
						energyTotalL = energyTotalR;
					} else if (!solidL && solidR) {	//left cell has a wall on the right
						densityR = densityL;
						//reflect along normal
						const velocityNL = normalX * velocityXL + normalY * velocityYL;
						velocityXR = velocityXL - 2 * normalX * velocityNL;
						velocityYR = velocityYL - 2 * normalY * velocityNL;
						energyTotalR = energyTotalL;
					}
				
					const energyKineticL = .5 * (velocityXL * velocityXL + velocityYL * velocityYL);
					const energyPotentialL = (xL - xmin) * _G.externalForceX + (yL - ymin) * _G.externalForceY;
					const energyThermalL = energyTotalL - energyKineticL - energyPotentialL;
					const pressureL = (this.gamma - 1) * densityL * energyThermalL;
					const speedOfSoundL = Math.sqrt(this.gamma * pressureL / densityL);
					const hTotalL = energyTotalL + pressureL / densityL;
					const roeWeightL = Math.sqrt(densityL);
				
					const energyKineticR = .5 * (velocityXR * velocityXR + velocityYR * velocityYR);
					const energyPotentialR = (xR - xmin) * _G.externalForceX + (yR - ymin) * _G.externalForceY;
					const energyThermalR = energyTotalR - energyKineticR - energyPotentialR;
					const pressureR = (this.gamma - 1) * densityR * energyThermalR;
					const speedOfSoundR = Math.sqrt(this.gamma * pressureR / densityR);
					const hTotalR = energyTotalR + pressureR / densityR;
					const roeWeightR = Math.sqrt(densityR);

					const denom = roeWeightL + roeWeightR;
					const velocityX = (roeWeightL * velocityXL + roeWeightR * velocityXR) / denom;
					const velocityY = (roeWeightL * velocityYL + roeWeightR * velocityYR) / denom;
					const hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) / denom;
					
					EulerEquationGodunovSolver.prototype.buildEigenstate(
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
	}
}

let eulerEquationSimulation = {
	methods : {
		'Burgers / Explicit' : EulerEquationBurgersExplicit.prototype,
		'Burgers / Backward Euler via Gauss Seidel' : EulerEquationBurgersBackwardEulerGaussSeidel.prototype,
		'Godunov / Explicit' : EulerEquationGodunovExplicit.prototype,
		'Roe / Explicit' : EulerEquationRoeExplicit.prototype,
	},
	initialConditions : {
		Sod : function() {
			let e = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * e];
					const y = this.x[1 + 2 * e];
					const rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
					let u = 0;
					let v = 0;
					if (useNoise) {
						u += (Math.random() - .5) * 2 * noiseAmplitude;
						v += (Math.random() - .5) * 2 * noiseAmplitude;
					}
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = 1;
					const energyTotal = energyKinetic + energyThermal + energyPotential;
					this.q[0 + 4 * e] = rho;
					this.q[1 + 4 * e] = rho * u; 
					this.q[2 + 4 * e] = rho * v; 
					this.q[3 + 4 * e] = rho * energyTotal; 
					this.solid[e] = 0;
					++e;
				}
			}
		},
		'Sod w/Cylinder' : function() {
			let e = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * e];
					const y = this.x[1 + 2 * e];
					const rho = ((x < (.7 * xmin + .3 * xmax) && y < (.7 * ymin + .3 * ymax)) ? 1 : .1);
					let u = 0;
					let v = 0;
					if (useNoise) {
						u += (Math.random() - .5) * 2 * noiseAmplitude;
						v += (Math.random() - .5) * 2 * noiseAmplitude;
					}
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = 1;
					const energyTotal = energyKinetic + energyThermal + energyPotential;
					this.q[0 + 4 * e] = rho;
					this.q[1 + 4 * e] = rho * u; 
					this.q[2 + 4 * e] = rho * v; 
					this.q[3 + 4 * e] = rho * energyTotal; 
			
					//insert a cylinder
					// ... with staggered rectangular boundaries
					const cx = .35 * xmin + .65 * xmax;
					const cy = .35 * ymin + .65 * ymax;
					const dx = x - cx;
					const dy = y - cy;
					const rSq = dx * dx + dy * dy;
					if (rSq < .1 * .1) {
						this.solid[e] = 1;
					} else {
						this.solid[e] = 0;
					}
					++e;
				}
			}
		},
		Wave : function() {
			const xmid = .5 * (xmin + xmax);
			const ymid = .5 * (ymin + ymax);
			const dg = .2 * (xmax - xmin);
			let e = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * e];
					const y = this.x[1 + 2 * e];
					const dx = x - xmid;
					const dy = y - ymid;
					const rho = 3 * Math.exp(-(dx * dx + dy * dy) / (dg * dg)) + .1;
					let u = 0; 
					let v = 0;
					if (useNoise) {
						u += (Math.random() - .5) * 2 * noiseAmplitude;
						v += (Math.random() - .5) * 2 * noiseAmplitude;
					}
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyThermal = 1;
					const energyTotal = energyKinetic + energyThermal + energyPotential;
					this.q[0 + 4 * e] = rho;
					this.q[1 + 4 * e] = rho * u;
					this.q[2 + 4 * e] = rho * v;
					this.q[3 + 4 * e] = rho * energyTotal;
					this.solid[e] = 0;
					++e;
				}
			}
		},
		//http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
		'Kelvin-Hemholtz' : function() {
			const xmid = .5 * (xmin + xmax);
			let e = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * e];
					const y = this.x[1 + 2 * e];
					const yInTheMiddle = y > (.75 * ymin + .25 * ymax) && y < (.25 * ymin + .75 * ymax);
					const rho = yInTheMiddle ? 2 : 1;
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
					const pressure = 2.5;
					const energyKinetic = .5 * (u * u + v * v);
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const energyTotal = pressure / ((this.gamma - 1) * rho) + energyKinetic + energyPotential;
					this.q[0 + 4 * e] = rho;
					this.q[1 + 4 * e] = rho * u; 
					this.q[2 + 4 * e] = rho * v;
					this.q[3 + 4 * e] = rho * energyTotal; 
					this.solid[e] = 0;
					++e;
				}
			}
			//TODO make it periodic on the left/right borders and reflecting on the top/bottom borders
		},
		//http://www.astro.virginia.edu/VITA/ATHENA/rt.html
		'Rayleigh-Taylor' : function() {
			const sizeX = xmax - xmin;
			const sizeY = ymax - ymin;
			const ymid = .5 * (ymin + ymax);
			let e = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					const x = this.x[0 + 2 * e];
					const y = this.x[1 + 2 * e];
					const yGreaterThanMid = y > ymid;
					const rho = yGreaterThanMid ? 2 : 1;
					let u = 0;
					const amplitude = .01;
					let v = .25 * amplitude * (1 + Math.cos(2 * Math.PI * x / sizeX)) * (1 + Math.cos(2 * Math.PI * y / sizeY));
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
					const energyPotential = (x - xmin) * _G.externalForceX + (y - ymin) * _G.externalForceY;
					const pressure = 2.5 - rho * energyPotential; 
					//let pressure = (this.gamma - 1) * rho * (2.5 - energyPotential);	//this fits with our non-external-force static ...
					/*
					pressure gradient zero ...
					dPressure / dy + rho * externalForce[side] = 0
					...integrate wrt y ...
					pressure + rho * externalForce.y * y = constant
					pressure = constant - rho * externalForce.y = constant - rho * energyPotential
					*/
					const energyKinetic = .5 * (u * u + v * v);
					/*
					energyTotal = energyKinetic + energyPotential + energyThermal
					energyTotal = .5 * (u*u + v*v) + energyPotential + pressure / (rho * (gamma - 1))
					*/
					const energyTotal = pressure / (rho * (this.gamma - 1)) + energyKinetic + energyPotential;
					this.q[0 + 4 * e] = rho;
					this.q[1 + 4 * e] = rho * u;
					this.q[2 + 4 * e] = rho * v;
					this.q[3 + 4 * e] = rho * energyTotal;
					this.solid[e] = 0;
					if (j == 0 || j == this.nx-1) this.solid[e] = 1;
					++e;
				}
			}
		}
	}
};

class HydroState { 
	constructor(args) {
		this.nx = args.size;
		this.cfl =.5;
		this.gamma = args.gamma;
		
		//x_i,j,dim: cell positions
		//0 <= i < this.nx
		this.x = new Float64Array(this.nx * this.nx * 2);
		let e = 0;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				this.x[e] = xmin + (xmax - xmin) * i / (this.nx-1); ++e;
				this.x[e] = ymin + (ymax - ymin) * j / (this.nx-1); ++e;
			}
		}
		
		//x_{i-1/2},{j-1/2},dim: interface positions
		//0 <= i,j < this.nx+1
		//0 <= dim < 2
		this.xi = new Float64Array((this.nx+1) * (this.nx+1) * 2);
		e = 0;
		for (let j = 0; j <= this.nx; ++j) {
			let J = Math.min(j, this.nx-1);
			for (let i = 0; i <= this.nx; ++i) {
				const I = Math.min(i, this.nx-1);
				if (i == 0) {
					this.xi[e] = 2 * this.x[0 + 2 * (1 + this.nx * J)] - this.x[0 + 2 * (2 + this.nx * J)];
				} else if (i == this.nx) {
					this.xi[e] = 2 * this.x[0 + 2 * (this.nx-1 + this.nx * J)] - this.x[0 + 2 * (this.nx-2 + this.nx * J)];
				} else {
					this.xi[e] = .5*(this.x[0 + 2 * (i + this.nx * J)] + this.x[0 + 2 * (i-1 + this.nx * J)]);
				}
				++e;
				if (j == 0) {
					this.xi[e] = 2 * this.x[1 + 2 * (I + this.nx * 1)] - this.x[1 + 2 * (I + this.nx * 2)];
				} else if (j == this.nx) {
					this.xi[e] = 2 * this.x[1 + 2 * (I + this.nx * (this.nx-1))] - this.x[1 + 2 * (I + this.nx * (this.nx-2))];
				} else {
					this.xi[e] = .5*(this.x[1 + 2 * (I + this.nx * j)] + this.x[1 + 2 * (I + this.nx * (j-1))]);
				}
				++e;
			}
		}

		//q_i,j,state: state vector, stored as q[state + 4 * (j + this.nx * i)]
		//q_i,j,0: density: rho
		//q_i,j,1: momentum: rho * u
		//q_i,j,2: momentum: rho * v
		//q_i,j,3: work: rho * e
		this.q = new Float64Array(this.nx * this.nx * 4);
		e = 0;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				for (let state = 0; state < 4; ++state) {
					this.q[e] = 0; ++e;
				}
			}
		}

		this.tmpq = [];
		for (let k = 0; k < 5; ++k) {
			this.tmpq[k] = new Float64Array(this.nx * this.nx * 4);
			e = 0;
			for (let j = 0; j < this.nx; ++j) {
				for (let i = 0; i < this.nx; ++i) {
					for (let state = 0; state < 4; ++state) {
						this.tmpq[k][e] = 0; ++e;
					}
				}
			}
		}


		//my first attempt at arbitrary boundaries ...
		this.solid = new Float64Array(this.nx * this.nx);
		e = 0;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				this.solid[e] = 0; ++e;
			}
		}

		eulerEquationSimulation.initialConditions.Sod.call(this);
		
		
		//TODO it is tempting to merge r, f, and ui into an edge structure
		//and associate them with the nodes on either side of them,
		//but then I would lose out on the 2nd-order contributions to the flux limiter.


		//f_{i-1/2},{j-1/2},side,state: cell flux
		this.flux = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4);
		e = 0;
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					for (let state = 0; state < 4; ++state) {
						this.flux[e] = 0; ++e;
					}
				}
			}
		}
	

		//used for Burgers
		
		
		//r_{i-1/2},{j-1/2},side,state	
		this.r = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4);;
		e = 0;
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					for (let state = 0; state < 4; ++state) {
						this.r[e] = 0; ++e;
					}
				}
			}
		}
		
		//only used with Burger's eqn advection code
		//u_{i-1/2},{j-1/2},dim: interface velocity
		this.ui = new Float64Array((this.nx+1) * (this.nx+1) * 2);
		e = 0;
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					this.ui[e] = 0; ++e;
				}
			}
		}

		//only used with Burgers pressure source code
		//p_i,j: pressure
		this.pressure = new Float64Array(this.nx * this.nx);
		e = 0;
		for (let j = 0; j < this.nx; ++j) {
			for (let i = 0; i < this.nx; ++i) {
				this.pressure[e] = 0; ++e;
			}
		}


		//used for Riemann
	

		//a_{i-1/2},{j-1/2},side,state,state
		this.interfaceMatrix = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4 * 4);
		this.interfaceEigenvalues = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4);
		this.interfaceEigenvectors = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4 * 4);
		this.interfaceEigenvectorsInverse = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4 * 4);
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					for (let stateJ = 0; stateJ < 4; ++stateJ) {
						for (let stateI = 0; stateI < 4; ++stateI) {
							//initialize to identity matrix
							this.interfaceMatrix[stateI + 4 * (stateJ + 4 * (side + 2 * (i + (this.nx+1) * j)))] = stateI == stateJ ? 1 : 0;
							this.interfaceEigenvectors[stateI + 4 * (stateJ + 4 * (side + 2 * (i + (this.nx+1) * j)))] = stateI == stateJ ? 1 : 0;
							this.interfaceEigenvectorsInverse[stateI + 4 * (stateJ + 4 * (side + 2 * (i + (this.nx+1) * j)))] = stateI == stateJ ? 1 : 0;
						}
						this.interfaceEigenvalues[stateJ + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
					}
				}
			}
		}

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

		//rTilde_{i-1/2},{j-1/2},side,state
		this.rTilde = new Float64Array((this.nx+1) * (this.nx+1) * 2 * 4);
		for (let j = 0; j <= this.nx; ++j) {
			for (let i = 0; i <= this.nx; ++i) {
				for (let side = 0; side < 2; ++side) {
					for (let state = 0; state < 4; ++state) {
						this.rTilde[state + 4 * (side + 2 * (i + (this.nx+1) * j))] = 0;
					}
				}
			}
		}

		//used for Backward Euler + Gauss-Seidel
		this.oldQ = new Float64Array(this.nx * this.nx * 4);

		//number of ghost cells
		this.nghost = 2;

		//solver configuration
		this.boundaryMethodTop = 'mirror';
		this.boundaryMethodLeft = 'mirror';
		this.boundaryMethodRight = 'mirror';
		this.boundaryMethodBottom = 'mirror';
		this.fluxMethod = 'superbee';
		this.algorithm = 'Roe / Explicit';
		this.explicitMethod = 'Euler';
		this.drawToScreenMethod = 'Density';
	}
	boundary() {
		boundaryMethods[this.boundaryMethodTop].top.call(this, this.nx, this.q);
		boundaryMethods[this.boundaryMethodLeft].left.call(this, this.nx, this.q);
		boundaryMethods[this.boundaryMethodRight].right.call(this, this.nx, this.q);
		boundaryMethods[this.boundaryMethodBottom].bottom.call(this, this.nx, this.q);
	}
	step(dt) {
		//apply boundary conditions
		this.boundary();
	
		//solve
		eulerEquationSimulation.methods[this.algorithm].step.call(this, dt);
	}
	update() {
		//do any pre-calcCFLTimestep preparation (Roe computes eigenvalues here)
		eulerEquationSimulation.methods[this.algorithm].initStep.call(this);
		
		//get timestep
		let dt;
		if (useCFL) {
			dt = eulerEquationSimulation.methods[this.algorithm].calcCFLTimestep.call(this);
		} else {
			dt = fixedDT;
		}
window.lastDT = dt;

		//do the update
		this.step(dt);
	}
}

class Hydro {
	constructor() {
		let size = Number(urlparams.get('size'));
		if (!size || !isFinite(size)) size = 200;
		let gamma = Number(urlparams.get('gamma'));
		if (!gamma || !isFinite(gamma)) gamma = 7/5;
		this.state = new HydroState({
			size : size,
			gamma : gamma
		});
	
		//geometry
		this.vertexPositions = new Float32Array(2*this.state.nx*this.state.nx);
		this.vertexStates = new Float32Array(this.state.nx*this.state.nx);
	
		//initialize geometry x and y coordinates once since they don't change
		const centerX = (xmax + xmin) / 2;
		const centerY = (ymax + ymin) / 2;
		let vIndex = 0;
		let xIndex = 0;
		const x = this.state.x;
		const nx = this.state.nx;
		for (let j = 0; j < nx; ++j) {
			for (let i = 0; i < nx; ++i) {
				this.vertexPositions[vIndex] = x[xIndex] - centerX; ++vIndex; ++xIndex;
				this.vertexPositions[vIndex] = x[xIndex] - centerY; ++vIndex; ++xIndex;
			}
		}
	}
	update() {
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
				ids.dataRangeFixedMin.value = this.lastDataMin;
				ids.dataRangeFixedMax.value = this.lastDataMax;
			}
		}
	}
}

const hydro = new Hydro();

function update() {
	//iterate
	hydro.update();
	waveVtxBuf.updateData(hydro.vertexPositions);
	waveStateBuf.updateData(hydro.vertexStates);
	//draw
	glutil.draw();
	requestAnimationFrame(update);
}

function onresize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	ids.content.style.height = (window.innerHeight - 50)+'px';
	glutil.resize();
}

function buildSelect(id, key, map) {
	let select = ids[id];
	for (let k in map) {
		let option = DOM('option', {text : k, appendTo : select});
		if (hydro.state[key] == k) {
			option.setAttribute('selected', 'true');
		}
	}
	select.addEventListener('change', e => {
		hydro.state[key] = select.value;
	});
	return select;
}

const sceneObjects = [];
const colorSchemes = {};

panel = ids.panel;

let panelContent = ids.content;
ids.menu.addEventListener('click', e => {
	if (hidden(panelContent)) {
		show(panelContent);
	} else {
		hide(panelContent);
	}
});

{
	let parent = ids.resetPanel;	//TODO specific panel for each simulation type
	let simulation = eulerEquationSimulation;
	for (let initialConditionName in simulation.initialConditions) {
		let method = simulation.initialConditions[initialConditionName];
		DOM('button', {
			text : 'Reset '+initialConditionName,
			click : function() {
				method.call(hydro.state);
			},
			appendTo : parent,
		});
		DOM('br', {appendTo : parent});
	}
}

ids.useNoise.addEventListener('change', e => {
	useNoise = ids.useNoise.checked;
});

{
	const select = ids.gridsize;
	[20, 50, 100, 200, 500, 1000].forEach(gridsize => {
		const option = DOM('option', {text : gridsize, appendTo : select});
		if (hydro.state.nx == gridsize) option.setAttribute('selected', 'true');
	});
	select.addEventListener('change', e => {
		const params = new URLSearchParams(urlparams);
		params.set('size', select.value);
		location.href = location.origin + location.pathname + '?' + params.toString();
	});
}

[
	'Top',
	'Left',
	'Right',
	'Bottom',
].forEach(sideName => {
	const selectName = 'boundaryMethod' + sideName;
	const constantName = 'boundary' + sideName + 'ConstantValue';
	const select = buildSelect(selectName, selectName, boundaryMethods);
	const onSelect = () => {
		if (hydro.state[selectName] == 'constant') {
			show(ids[constantName]);
		} else {
			hide(ids[constantName]);
		}
	}
	select.addEventListener('change', onSelect);
	onSelect();
});

buildSelect('fluxMethod', 'fluxMethod', fluxMethods);
buildSelect('drawToScreenMethod', 'drawToScreenMethod', drawToScreenMethods);
buildSelect('Euler_algorithm', 'algorithm', eulerEquationSimulation.methods);

buildSelect('explicitMethod', 'explicitMethod', explicitMethods);

[
	'externalForceX',
	'externalForceY',
	'boundaryTopConstantValue',
	'boundaryLeftConstantValue',
	'boundaryRightConstantValue',
	'boundaryBottomConstantValue',
].forEach(varName => {
	const thiz = ids[varName];
	thiz.value = _G[varName];
	thiz.addEventListener('change', e => {
		const v = Number(thiz.value);
		if (!isFinite(v)) return;	//NaN
		_G[varName] = v;
	});
});

hydro.lastDataMin = Number(ids.dataRangeFixedMin.value);
hydro.lastDataMax = Number(ids.dataRangeFixedMax.value);
hydro.updateLastDataRange = true;
ids.dataRangeScaleNormalized.addEventListener('change', e => {
	if (!ids.dataRangeScaleNormalized.checked) return;
	hydro.updateLastDataRange = true;
});
ids.dataRangeScaleFixed.addEventListener('change', e => {
	if (!ids.dataRangeScaleFixed.checked) return;
	hydro.updateLastDataRange = false;
	hydro.lastDataMin = Number(ids.dataRangeFixedMin.value);
	hydro.lastDataMax = Number(ids.dataRangeFixedMax.value);
});
ids.dataRangeFixedMin.addEventListener('change', e => {
	if (hydro.updateLastDataRange) return;
	hydro.lastDataMin = Number(ids.dataRangeFixedMin.value);
});
ids.dataRangeFixedMax.addEventListener('change', e => {
	if (hydro.updateLastDataRange) return;
	hydro.lastDataMax = Number(ids.dataRangeFixedMax.value);
});

ids.timeStepCFLBased.addEventListener('change', e => {
	if (!ids.timeStepCFLBased.checked) return;
	useCFL = true;
});
ids.timeStepCFL.value = hydro.state.cfl;
ids.timeStepCFL.addEventListener('change', e => {
	let v = Number(ids.timeStepCFL.value);
	if (!isFinite(v)) return;
	hydro.state.cfl = v;
});
ids.timeStepFixed.addEventListener('change', e => {
	if (!ids.timeStepFixed.checked) return;
	useCFL = false;
});
ids.timeStepValue.value = fixedDT;
ids.timeStepValue.addEventListener('change', e => {
	let v = Number(ids.timeStepValue.value);
	if (!isFinite(v)) return;	//stupid javascript ... convert anything it doesn't understand to NaNs...
	fixedDT = v;
});

canvas = DOM('canvas', {
	css : {
		left : 0,
		top : 0,
		position : 'absolute',
		userSelect : 'none',
	},
	prependTo : document.body,
});

try {
	glutil = new GLUtil({canvas:canvas});
	gl = glutil.context;
} catch (e) {
	removeFromParent(panel);
	removeFromParent(canvas);
	show(ids.webglfail);
	throw e;
}
glutil.import('Gradient', makeGradient);

gl.enable(gl.DEPTH_TEST);
glutil.view.ortho = true;
glutil.view.zNear = -1;
glutil.view.zFar = 1;
glutil.view.fovY = 125 / 200 * (xmax - xmin);

colorSchemes.Heat = new glutil.Gradient.GradientTexture({
	width : 256, 
	colors : [
		[0, 0, 0],
		[0, 0, 1],
		[1, 1, 0],
		[1, 0, 0],
	],
	dontRepeat : true
});

let isobarSize = 16;
let isobarData = new Uint8Array(isobarSize);
for (let i = 1; i < isobarSize; i += 2) {
	isobarData[i] = 255;
}
colorSchemes['B&W'] = new glutil.Texture2D({
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

for (let k in colorSchemes) {
	DOM('option', {
		value : k,
		text : k,
		appendTo : ids.colorScheme,
	});
}
ids.colorScheme.addEventListener('change', e => {
	const k = ids.colorScheme.value;
	colorSchemes[k].bind();
	//const v = colorSchemes[k];
	//sceneObjects.forEach(sceneObject => {
	//gl.bindTexture(gl.TEXTURE_2D, v.obj);
	//});
});

const shader = new glutil.Program({
	vertexCode : `
in vec2 vertex;
in float state;
out float statev;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	statev = state;
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
}
`,
	fragmentCode : `
in float statev;
uniform sampler2D tex;
out vec4 fragColor;
void main() {
	fragColor = texture(tex, vec2(statev, .5)); 
}
`,
});

//make grid
hydro.update();
waveVtxBuf = new glutil.ArrayBuffer({
	dim : 2,
	data : hydro.vertexPositions,
	usage : gl.DYNAMIC_DRAW
});
waveStateBuf = new glutil.ArrayBuffer({
	dim : 1,
	data : hydro.vertexStates,
	usage : gl.DYNAMIC_DRAW
});
for (let j = 0; j < hydro.state.nx-1; ++j) {
	let indexes = [];
	for (let i = 0; i < hydro.state.nx; ++i) {
		indexes.push(i + j*hydro.state.nx);
		indexes.push(i + (j+1)*hydro.state.nx);
	}
	sceneObjects.push(new glutil.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		indexes : new glutil.ElementArrayBuffer({
			data : new Uint32Array(indexes)
		}),
		attrs : {
			vertex : waveVtxBuf,
			state : waveStateBuf
		},
		uniforms : {
			tex : 0
		},
		shader : shader
	}));
}

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
		glutil.view.pos[0] -= dx / canvas.width * 2 * (aspectRatio * glutil.view.fovY);
		glutil.view.pos[1] += dy / canvas.height * 2 * glutil.view.fovY;
		glutil.updateProjection();
	},
	zoom : function(zoomChange) {
		dragging = true;
		let scale = Math.exp(-zoomFactor * zoomChange);
		glutil.view.fovY *= scale 
		glutil.updateProjection();
	}
});

//start it off
window.addEventListener('resize', onresize);
onresize();
update();
