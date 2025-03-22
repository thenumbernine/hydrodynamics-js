/*
1D ADM
2D: Burger at least
2D ADM
2D unstructured
*/

import {Br, Button, Canvas, Option} from '/js/dom.js';
import {assert} from '/js/util.js';
import {svd} from './svd.js';
import {getIDs, removeFromParent, show, hide, hidden} from '/js/util.js';
import {GLUtil} from '/js/gl-util.js';
import {Mouse3D} from '/js/mouse3d.js';
import {makeExplicitMethods} from './explicit_methods.js';
import {fluxMethods} from './flux_methods.js';

const ids = getIDs();
window.ids = ids;

const urlparams = new URLSearchParams(location.search);

let gl;
let glutil;
let mouse;
let panel;
let canvas;
let xmin = -1;
let xmax = 1; 
let ymin = -1;
let ymax = 2;
let gridstep = .1;
let useCFL = true;
let fixedDT = .2;
let pause = false;
let gaussSeidelIterations = 20;
let plainShader;
let graphShader;
let gridObj;

function mat33invert(out, a) {
	let det = a[0][0] * a[1][1] * a[2][2]
			+ a[1][0] * a[2][1] * a[0][2]
			+ a[2][0] * a[0][1] * a[1][2]
			- a[2][0] * a[1][1] * a[0][2]
			- a[1][0] * a[0][1] * a[2][2]
			- a[0][0] * a[2][1] * a[1][2];
	if (det == 0) {
		console.log('a', a);
		console.log('singular!');
assert(false); //get a proper stack trace since throw doesn't
		return;
	}
	let invDet = 1 / det;
	for (let j = 0; j < 3; ++j) {
		let j1 = (j + 1) % 3;
		let j2 = (j + 2) % 3;
		for (let i = 0; i < 3; ++i) {
			let i1 = (i + 1) % 3;
			let i2 = (i + 2) % 3;
			out[i][j] = invDet * (a[j1][i1] * a[j2][i2] - a[j1][i2] * a[j2][i1]);
		}
	}
}

let copyState = function(srcQ, destQ) {
	for (let i = 0; i < srcQ.length; ++i) {
		for (let j = 0; j < 3; ++j) {
			destQ[i][j] = srcQ[i][j];
		}
	}
	return destQ;
};

let addMulState = function(to, from, scalar) {
	for (let i = 0; i < to.length; ++i) {
		for (let j = 0; j < 3; ++j) {
			to[i][j] += scalar * from[i][j];
		}
	}
};

const explicitMethods = makeExplicitMethods({
	copyState : copyState,
	addMulState : addMulState,
});

let boundaryMethods = {
	periodic : function(nx,q) {
		q[0][0] = q[nx-4][0];
		q[1][0] = q[nx-3][0];
		q[nx-2][0] = q[2][0];
		q[nx-1][0] = q[3][0];
		q[0][1] = q[nx-4][1];
		q[1][1] = q[nx-3][1];
		q[nx-2][1] = q[2][1];
		q[nx-1][1] = q[3][1];
		q[0][2] = q[nx-4][2];
		q[1][2] = q[nx-3][2];
		q[nx-2][2] = q[2][2];
		q[nx-1][2] = q[3][2];
	},
	mirror : function(nx,q) {
		q[0][0] = q[3][0];
		q[1][0] = q[2][0];
		q[nx-2][0] = q[nx-3][0];
		q[nx-1][0] = q[nx-4][0];
		q[0][1] = -q[3][1];
		q[1][1] = -q[2][1];
		q[nx-2][1] = -q[nx-3][1];
		q[nx-1][1] = -q[nx-4][1];
		q[0][2] = q[3][2];
		q[1][2] = q[2][2];
		q[nx-2][2] = q[nx-3][2];
		q[nx-1][2] = q[nx-4][2];
	},
	/*constant : function(nx,q) {
		q[0][0] = 0;
		q[1][0] = 0;
		q[nx-2][0] = 0;
		q[nx-1][0] = 0;
		q[0][1] = 0;
		q[1][1] = 0;
		q[nx-2][1] = 0;
		q[nx-1][1] = 0;
		q[0][2] = 0;
		q[1][2] = 0;
		q[nx-2][2] = 0;
		q[nx-1][2] = 0;
	},
	*/
	freeflow : function(nx,q) {
		q[0][0] = q[1][0] = q[2][0];
		q[nx-1][0] = q[nx-2][0] = q[nx-3][0];
		q[0][1] = q[1][1] = q[2][1];
		q[nx-1][1] = q[nx-2][1] = q[nx-3][1];
		q[0][2] = q[1][2] = q[2][2];
		q[nx-1][2] = q[nx-2][2] = q[nx-3][2];
	}
};

function eulerGetPrimitives(i) {
	let rho = this.q[i][0];
	let u = this.q[i][1] / rho;
	let e = this.q[i][2] / rho;
	return [rho, u, e];
}

class EulerEquationBurgersSolver {
	initStep() {}
	calcCFLTimestep() {
		let mindum = undefined;
		for (let i = 0; i < this.nx; ++i) {
			let u = this.q[i][1] / this.q[i][0];
			let energyTotal = this.q[i][2] / this.q[i][0];
			let energyKinetic = .5 * u * u;
			let energyInternal = energyTotal - energyKinetic;
			let speedOfSound = Math.sqrt(this.gamma * (this.gamma - 1) * energyInternal);
			let dx = this.xi[i+1] - this.xi[i];
			let dum = dx / (speedOfSound + Math.abs(u));
			if (mindum === undefined || dum < mindum) mindum = dum;
		}
		//if (mindum != mindum) throw 'nan';
		return this.cfl * mindum;
	}
	getPrimitives(...args) {
		return eulerGetPrimitives.call(this, ...args);
	}
}

/*
time derivative + advection component + source term = 0
d/dt (rho) + d/dx(rho u) = 0
d/dt (rho u) + d/dx(rho u&2) + d/dx (P) = 0
d/dt (rho e_total) + d/dx (rho e_total u) + d/dx (P u) = 0
*/
class EulerEquationBurgersExplicit extends EulerEquationBurgersSolver {
	step(dt) {
		//boundary precedes this call
		explicitMethods[this.explicitMethod].call(this, dt, EulerEquationBurgersExplicit.prototype.integrateFlux);
		
		//boundary again
		this.boundary();
		
		explicitMethods[this.explicitMethod].call(this, dt, EulerEquationBurgersExplicit.prototype.integrateMomentumDiffusion);
		explicitMethods[this.explicitMethod].call(this, dt, EulerEquationBurgersExplicit.prototype.integrateWorkDiffusion);
	}
	integrateFlux(dt, dq_dt) {
		assert(this.x.length == this.nx);
		assert(this.xi.length == this.nx + 1);
		assert(this.q.length == this.nx);
		assert(this.ui.length == this.nx + 1);
	
		//get velocity at interfaces from state
		for (let ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
			this.ui[ix] = .5 * (this.q[ix][1] / this.q[ix][0] + this.q[ix-1][1] / this.q[ix-1][0]);
		}
		this.ui[0] = this.ui[this.nx] = 0;

		//compute flux and advect for each state vector
		for (let j = 0; j < 3; ++j) {
			//r_{i-1/2} flux limiter
			for (let i = this.nghost; i < this.nx+1-this.nghost; ++i) {
				let dq = this.q[i][j] - this.q[i-1][j];
				if (Math.abs(dq) > 0) {
					if (this.ui[i] >= 0) {
						this.r[i][j] = (this.q[i-1][j] - this.q[i-2][j]) / dq;
					} else {
						this.r[i][j] = (this.q[i+1][j] - this.q[i][j]) / dq;
					}
				} else {
					this.r[i][j] = 0;
				}
			}
			this.r[0][j] = this.r[1][j] = this.r[this.nx-1][j] = this.r[this.nx][j] = 0;

			//construct flux:
			for (let i = this.nghost-1; i < this.nx+this.nghost-2; ++i) {
				//flux limiter
				let phi = fluxMethods[this.fluxMethod](this.r[i][j]);
				if (this.ui[i] >= 0) {
					this.interfaceFlux[i][j] = this.ui[i] * this.q[i-1][j];
				} else {
					this.interfaceFlux[i][j] = this.ui[i] * this.q[i][j];
				}
				let delta = phi * (this.q[i][j] - this.q[i-1][j]);
				let dx = this.x[i] - this.x[i-1];
				this.interfaceFlux[i][j] += delta * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
			}
			this.interfaceFlux[0][j] = this.interfaceFlux[this.nx][j] = 0;

			//update cells
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				dq_dt[i][j] = -(this.interfaceFlux[i+1][j] - this.interfaceFlux[i][j]) / (this.xi[i+1] - this.xi[i]);
			}
			for (let i = 0; i < this.nghost; ++i) {
				dq_dt[i] = [0,0,0];
				dq_dt[this.nx-i-1] = [0,0,0];
			}
		}			
	}
	//compute pressure
	integrateMomentumDiffusion(dt, dq_dt) {
		for (let i = 0; i < this.nx; ++i) {
			let u = this.q[i][1] / this.q[i][0];
			let energyTotal = this.q[i][2] / this.q[i][0];
			let energyKinetic = .5 * u * u;
			let energyInternal = energyTotal - energyKinetic;
			this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyInternal;
		}

		//apply momentum diffusion = pressure
		for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
			dq_dt[i][0] = 0;
			dq_dt[i][1] = -(this.pressure[i+1] - this.pressure[i-1]) / (this.x[i+1] - this.x[i-1]);
			dq_dt[i][2] = 0;
		}
	}
	integrateWorkDiffusion(dt, dq_dt) {
		//apply work diffusion = momentum
		for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
			let u_inext = this.q[i+1][1] / this.q[i+1][0];
			let u_iprev = this.q[i-1][1] / this.q[i-1][0];
			dq_dt[i][0] = 0;
			dq_dt[i][1] = 0;
			dq_dt[i][2] = -(this.pressure[i+1] * u_inext - this.pressure[i-1] * u_iprev) / (this.x[i+1] - this.x[i-1]);
		}
	}
}

/*
TODO 3x3 block tridiagonal thomas algorithm
until then, Gauss-Seidel

Linearizing our relationship between current and next timesteps.
Treating the flux limiter and the interface velocity as constants when I could consider them in terms of q's. 
I get timesteps of .7, when .3 or so is what the max CFL timestep for explicit was giving me,
 but still see a lot more oscillations in the system.
*/
class EulerEquationBurgersBackwardEulerGaussSeidel extends EulerEquationBurgersSolver {
	step(dt) {
		for (let i = 0; i < this.nx; ++i) {
			for (let j = 0; j < 3; ++j) {
				this.oldQ[i][j] = this.q[i][j];
			}
		}
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			//get velocity at interfaces from state
			for (let ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
				this.ui[ix] = .5 * (this.q[ix][1] / this.q[ix][0] + this.q[ix-1][1] / this.q[ix-1][0]);
			}
			this.ui[0] = this.ui[this.nx] = 0;

			//compute flux and advect for each state vector
			for (let j = 0; j < 3; ++j) {
				//r_{i-1/2} flux limiter
				for (let i = this.nghost; i < this.nx+this.nghost-3; ++i) {
					let dq = this.q[i][j] - this.q[i-1][j];
					if (Math.abs(dq) > 0) {
						if (this.ui[i] >= 0) {
							this.r[i][j] = (this.q[i-1][j] - this.q[i-2][j]) / dq;
						} else {
							this.r[i][j] = (this.q[i+1][j] - this.q[i][j]) / dq;
						}
					} else {
						this.r[i][j] = 0;
					}
				}
				this.r[0][j] = this.r[1][j] = this.r[this.nx-1][j] = this.r[this.nx][j] = 0;

				//update cells
				
				for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
					let dx = this.xi[i+1] - this.xi[i];
					
					let dflux_qip = 0;
					let dflux_qi = 0;
					let dflux_qin = 0;
					
					//flux limiter
					let phi = fluxMethods[this.fluxMethod](this.r[i][j]);
					if (this.ui[i] >= 0) {
						dflux_qip -= this.ui[i];
					} else {
						dflux_qi -= this.ui[i];
					}
					dflux_qi -= phi * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
					dflux_qip += phi * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
					
					phi = fluxMethods[this.fluxMethod](this.r[i+1][j]);
					if (this.ui[i+1] >= 0) {
						dflux_qi += this.ui[i+1];
					} else {
						dflux_qin += this.ui[i+1];
					}
					dflux_qin += phi * .5 * Math.abs(this.ui[i+1]) * (1 - Math.abs(this.ui[i+1] * dt / dx));
					dflux_qi -= phi * .5 * Math.abs(this.ui[i+1]) * (1 - Math.abs(this.ui[i+1] * dt / dx));
					
					this.q[i][j] = (this.oldQ[i][j] - dt / dx * (
						dflux_qip * this.q[i-1][j] 
						+ dflux_qin * this.q[i+1][j])) / (1 + dt / dx * dflux_qi);
				}
			}
		}	
	
		// only needed for Burgers	
		
		//boundary again
		this.boundary();

		//apply momentum diffusion = pressure
		for (let i = 0; i < this.nx; ++i) {
			for (let j = 0; j < 3; ++j) {
				this.oldQ[i][j] = this.q[i][j];
			}
		}
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			for (let i = 0; i < this.nx; ++i) {
				let u = this.q[i][1] / this.q[i][0];
				let energyTotal = this.q[i][2] / this.q[i][0];
				let energyKinetic = .5 * u * u;
				let energyInternal = energyTotal - energyKinetic;
				this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyInternal;
			}
			
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				this.q[i][1] = this.oldQ[i][1] - dt * (this.pressure[i+1] - this.pressure[i-1]) / (this.x[i+1] - this.x[i-1]);
			}
		}

		//apply work diffusion = momentum
		for (let i = 0; i < this.nx; ++i) {
			for (let j = 0; j < 3; ++j) {
				this.oldQ[i][j] = this.q[i][j];
			}
		}
		for (let iter = 0; iter < gaussSeidelIterations; ++iter) {
			for (let i = 0; i < this.nx; ++i) {
				let u = this.q[i][1] / this.q[i][0];
				let energyTotal = this.q[i][2] / this.q[i][0];
				let energyKinetic = .5 * u * u;
				let energyInternal = energyTotal - energyKinetic;
				this.pressure[i] = (this.gamma - 1) * this.q[i][0] * energyInternal;
			}
			
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				let u_inext = this.q[i+1][1] / this.q[i+1][0];
				let u_iprev = this.q[i-1][1] / this.q[i-1][0];
				this.q[i][2] = this.oldQ[i][2] - dt * (this.pressure[i+1] * u_inext - this.pressure[i-1] * u_iprev) / (this.x[i+1] - this.x[i-1]);
			}
		}
	}
}

class EulerEquationBurgersBackwardEulerTridiagonal extends EulerEquationBurgersSolver {
	step(dt) {
		//get velocity at interfaces from state
		for (let ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
			this.ui[ix] = .5 * (this.q[ix][1] / this.q[ix][0] + this.q[ix-1][1] / this.q[ix-1][0]);
		}
		this.ui[0] = this.ui[this.nx] = 0;

		//compute flux and advect for each state vector
		for (let j = 0; j < 3; ++j) {
			//r_{i-1/2} flux limiter
			for (let i = this.nghost; i < this.nx+this.nghost-3; ++i) {
				let dq = this.q[i][j] - this.q[i-1][j];
				if (Math.abs(dq) > 0) {
					if (this.ui[i] >= 0) {
						this.r[i][j] = (this.q[i-1][j] - this.q[i-2][j]) / dq;
					} else {
						this.r[i][j] = (this.q[i+1][j] - this.q[i][j]) / dq;
					}
				} else {
					this.r[i][j] = 0;
				}
			}
			this.r[0][j] = this.r[1][j] = this.r[this.nx-1][j] = this.r[this.nx][j] = 0;

			//update cells
			
			for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
				let dx = this.xi[i+1] - this.xi[i];
				
				let dflux_qip = 0;
				let dflux_qi = 0;
				let dflux_qin = 0;
				
				//flux limiter
				let phi = fluxMethods[this.fluxMethod](this.r[i][j]);
				if (this.ui[i] >= 0) {
					dflux_qip -= this.ui[i];
				} else {
					dflux_qi -= this.ui[i];
				}
				dflux_qi -= phi * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
				dflux_qip += phi * .5 * Math.abs(this.ui[i]) * (1 - Math.abs(this.ui[i] * dt / dx));
				
				phi = fluxMethods[this.fluxMethod](this.r[i+1][j]);
				if (this.ui[i+1] >= 0) {
					dflux_qi += this.ui[i+1];
				} else {
					dflux_qin += this.ui[i+1];
				}
				dflux_qin += phi * .5 * Math.abs(this.ui[i+1]) * (1 - Math.abs(this.ui[i+1] * dt / dx));
				dflux_qi -= phi * .5 * Math.abs(this.ui[i+1]) * (1 - Math.abs(this.ui[i+1] * dt / dx));
				
				this.q[i][j] = (this.q[i][j] - dt / dx * (
					dflux_qip * this.q[i-1][j] 
					+ dflux_qin * this.q[i+1][j])) / (1 + dt / dx * dflux_qi);
			}
		}

		// only needed for Burgers	
		
		//boundary again
		this.boundary();

		let a = [];	//lower band
		let b = [];	//diagonal
		let c = [];	//upper band

		//apply momentum diffusion = pressure
		for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
			let dtOverTwoDx = dt / (this.x[i+1] - this.x[i-1]);
			this.q[i][1] = this.q[i][1] 
				- dtOverTwoDx * (this.gamma - 1) * (this.q[i+1][2] - .5 * this.q[i+1][1] * this.q[i+1][1] / this.q[i+1][0]) 
				+ dtOverTwoDx * (this.gamma - 1) * (this.q[i-1][2] - .5 * this.q[i-1][1] * this.q[i-1][1] / this.q[i-1][0]);
		}

		//apply work diffusion = momentum
		for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
			let dtOverTwoDx = dt / (this.x[i+1] - this.x[i-1]);
			this.q[i][2] = this.q[i][2] 
				- dtOverTwoDx * this.pressure[i+1] * this.q[i+1][1] / this.q[i+1][0]
				+ dtOverTwoDx * this.pressure[i-1] * this.q[i-1][1] / this.q[i-1][0];
		}
	}
}

class GodunovSolver {
	initStep() {
		let nx = this.nx;
		/* linear */
		for (let i = 0; i < nx; ++i) {
			for (let j = 0; j < 3; ++j) {
				this.qL[i][j] = this.q[i][j];
				this.qR[i][j] = this.q[i][j];
			}
		}
		for (let ix = 1; ix < nx; ++ix) {
			for (let j = 0; j < 3; ++j) {
				let iqL = this.qR[ix-1][j];
				let iqR = this.qL[ix][j];
				this.qMid[ix][j] = (iqL + iqR) * .5;
			}
		}
		/**/
		/* PPM * /
		for (let j = 0; j < 3; ++j) {
			this.qL[0][j] = this.q[0][j];
			this.qR[0][j] = this.q[0][j];
			this.qL[1][j] = this.q[1][j];
			this.qR[1][j] = this.q[1][j];
			this.qL[nx-2][j] = this.q[nx-2][j];
			this.qR[nx-2][j] = this.q[nx-2][j];
			this.qL[nx-1][j] = this.q[nx-1][j];
			this.qR[nx-1][j] = this.q[nx-1][j];
		}
		for (let ix = 2; ix < nx-1; ++ix) {
			for (let j = 0; j < 3; ++j) {
				this.qMid[ix][j] = 6/12 * (this.q[ix-1][j] + this.q[ix][j])
								- 1/12 * (this.q[ix-2][j] + this.q[ix+1][j]);
				this.qR[ix-1][j] = this.qMid[ix][j];
				this.qL[ix][j] = this.qMid[ix][j];
			}
		}	
		for (let j = 0; j < 3; ++j) {
			this.qMid[1][j] = (this.qR[0][j] + this.qL[1][j]) * .5;
			this.qMid[nx-1][j] = (this.qR[nx-2][j] + this.qL[nx-1][j]) * .5;
		}
		for (let i = 0; i < nx; ++i) {
			let q = this.q[i];
			let qR = this.qR[i];
			let qL = this.qL[i];
			for (let j = 0; j < 3; ++j) {
				if ((qR[j] - q[j]) * (q[j] - qL[j]) <= 0) {
					qL[j] = q[j];
					qR[j] = q[j];
				} else {
					let a = (qR[j] - qL[j]) * (q[j] - .5 * (qL[j] + qR[j]));
					let b = (qR[j] - qL[j]) * (qR[j] - qL[j]) / 6;
					if (a > b) {
						qL[j] = 3 * q[j] - 2 * qR[j];
					}
					if (a < -b) {
						qR[j] = 3 * q[j] - 2 * qL[j];
					}
				}
			}
		}
		/**/
	}
	//update first calls initStep and then aclls calcCFLTimestep
	//initStep is abstract.  it is where the eigen decomposition takes place.
	/*
	store eigenvalues and eigenvectors of interfaces
	use the lambdas to calc the DT based on CFL
	*/
	calcCFLTimestep() {

		let mindum = undefined;
		for (let i = 1; i < this.nx; ++i) {
			let maxLambda = Math.max(0, this.interfaceEigenvalues[i][0], this.interfaceEigenvalues[i][1], this.interfaceEigenvalues[i][2]);
			let minLambda = Math.min(0, this.interfaceEigenvalues[i+1][0], this.interfaceEigenvalues[i+1][1], this.interfaceEigenvalues[i+1][2]);
			let dum = (this.xi[i+1] - this.xi[i]) / (maxLambda - minLambda);
			if (mindum === undefined || dum < mindum) mindum = dum;
		}
		//if (mindum != mindum) throw 'nan';
		return this.cfl * mindum;
	}
	step(dt) {
		let deriv = GodunovSolver.prototype.calcDerivative;
		explicitMethods[this.explicitMethod].call(this, dt, deriv);
	}
	calcDerivative(dt, dq_dt) {
		for (let ix = 1; ix < this.nx; ++ix) {
			let iqL = this.qR[ix-1];
			let iqR = this.qL[ix];
			for (let j = 0; j < 3; ++j) {
				//the change in state represented in interface eigenbasis
				this.interfaceDeltaQTilde[ix][j] = 
					this.interfaceEigenvectorsInverse[ix][0][j] * (iqR[0] - iqL[0])
					+ this.interfaceEigenvectorsInverse[ix][1][j] * (iqR[1] - iqL[1])
					+ this.interfaceEigenvectorsInverse[ix][2][j] * (iqR[2] - iqL[2]);
			}
		}
		
		for (let j = 0; j < 3; ++j) {
			this.interfaceDeltaQTilde[0][j] = 0;
			this.interfaceDeltaQTilde[this.nx][j] = 0;
		}
		
		for (let ix = this.nghost; ix < this.nx-this.nghost+1; ++ix) {
			for (let j = 0; j < 3; ++j) {
				let interfaceDeltaQTilde = this.interfaceDeltaQTilde[ix][j];
				if (Math.abs(interfaceDeltaQTilde) > 0) {
					if (this.interfaceEigenvalues[ix][j] >= 0) {
						this.rTilde[ix][j] = this.interfaceDeltaQTilde[ix-1][j] / interfaceDeltaQTilde;
					} else {
						this.rTilde[ix][j] = this.interfaceDeltaQTilde[ix+1][j] / interfaceDeltaQTilde;
					}
				} else {
					this.rTilde[ix][j] = 0;
				}
			}
		}

		//..and keep the boundary r's zero	
		for (let j = 0; j < 3; ++j) {
			this.rTilde[0][j] = this.rTilde[1][j] = this.rTilde[this.nx-1][j] = this.rTilde[this.nx][j] = 0;
		}
		/*
		for (let ix = 0; ix < this.nghost; ++ix) {
			for (let j = 0; j < 3; ++j) {
				this.rTilde[ix][j] = 0;
				this.rTilde[this.nx-ix][j] = 0;
			}	
		}
		*/


		//transform cell q's into cell qTilde's (eigenspace)
		// ... so q_{i-1/2}L = q_{i-1}, q_{i-1/2}R = q_i
		// qTilde_{i-1/2}L = E_{i-1/2}^-1 q_{i-1}, qTilde_{i-1/2}R = E_{i-1/2}^-1 q_i
		//use them to detemine qTilde's at boundaries
		//use them (and eigenvalues at boundaries) to determine fTilde's at boundaries
		//use them (and eigenvectors at boundaries) to determine f's at boundaries
		//use them to advect, like good old fluxes advect
		let fluxTilde = [];
		let fluxAvg = [];

		//qi[ix] = q_{i-1/2} lies between q_{i-1} = q[i-1] and q_i = q[i]
		//(i.e. qi[ix] is between q[ix-1] and q[ix])
		//Looks good according to "Riemann Solvers and Numerical Methods for Fluid Dynamics," Toro, p.191
		for (let ix = this.nghost-1; ix < this.nx+this.nghost-2; ++ix) {
			//simplification: rather than E * L * E^-1 * q, just do A * q for A the original matrix
			//...and use that on the flux L & R avg (which doesn't get scaled in eigenvector basis space
			
			//if I wasn't doing all the above slope-limited stuff, this would suffice:
			//however if I'm not using this then I don't need to store the interface jacobian matrix
			for (let j = 0; j < 3; ++j) {
				fluxAvg[j] = 
					this.interfaceMatrix[ix][0][j] * this.qMid[ix][0]
					+ this.interfaceMatrix[ix][1][j] * this.qMid[ix][1]
					+ this.interfaceMatrix[ix][2][j] * this.qMid[ix][2];
			}

			//calculate flux
			for (let j = 0; j < 3; ++j) {
				let theta = 0;
				if (this.interfaceEigenvalues[ix][j] >= 0) {
					theta = 1;
				} else {
					theta = -1;
				}
				
				let phiTilde = fluxMethods[this.fluxMethod](this.rTilde[ix][j]);
				let dx = this.xi[ix] - this.xi[ix-1];
				let epsilon = this.interfaceEigenvalues[ix][j] * dt / dx;
				
				//interfaceFlux[ix][k] = fluxTilde[ix][j] * interfaceEigenvectors[ix][k][j]
				//flux in eigenvector basis is the q vector transformed by the inverse then scaled by the eigenvalue
				//should the eigenvalue be incorperated here, after flux limiter is taken into account, or beforehand?
				//1D says after, but notes say before ...
				let deltaFluxTilde = this.interfaceEigenvalues[ix][j] * this.interfaceDeltaQTilde[ix][j];
				
				fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phiTilde * (epsilon - theta));
			}

			//reproject fluxTilde back into q
			for (let j = 0; j < 3; ++j) {
				this.interfaceFlux[ix][j] = fluxAvg[j]
					+ this.interfaceEigenvectors[ix][0][j] * fluxTilde[0] 
					+ this.interfaceEigenvectors[ix][1][j] * fluxTilde[1] 
					+ this.interfaceEigenvectors[ix][2][j] * fluxTilde[2];
			}
		}
		
		//zero boundary flux
		for (let j = 0; j < 3; ++j) {
			this.interfaceFlux[0][j] = this.interfaceFlux[this.nx][j] = 0;
		}

		//update cells
		for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
			for (let j = 0; j < 3; ++j) {
				dq_dt[i][j] = -(this.interfaceFlux[i+1][j] - this.interfaceFlux[i][j]) / (this.xi[i+1] - this.xi[i]);
			}
		}
		for (let i = 0; i < this.nghost; ++i) {
			dq_dt[i] = [0,0,0];
			dq_dt[this.nx-i-1] = [0,0,0];
		}
	}
}
			
window.maximalError = 0;
window.maximalEigenError = 0;

class EulerEquationGodunovSolver extends GodunovSolver {}
EulerEquationGodunovSolver.prototype.buildEigenstate = {
	Analytic : {
		calcFlux : function(matrix, velocity, hTotal, gamma, eTotal) {
			//flux jacobian matrix, listed per column
			matrix[0][0] = 0;
			matrix[0][1] = (gamma - 3) / 2 * velocity * velocity;
			matrix[0][2] = -velocity * (gamma * eTotal + (1 - gamma) * velocity * velocity); 
			matrix[1][0] = 1;
			matrix[1][1] = (3 - gamma) * velocity;
			matrix[1][2] = hTotal + (1 - gamma) * velocity * velocity;
			matrix[2][0] = 0;
			matrix[2][1] = gamma - 1;
			matrix[2][2] = gamma * velocity;
		},
		calcEigenvalues : function(eigenvalues, velocity, speedOfSound) {
			//eigenvalues: min, mid, max
			eigenvalues[0] = velocity - speedOfSound;
			eigenvalues[1] = velocity;
			eigenvalues[2] = velocity + speedOfSound;
		},
		calcAll : function(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, hTotal, gamma, eTotal) {
			//calculate matrix & eigenvalues & vectors at interface from state at interface
			let speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * velocity * velocity));	

			EulerEquationGodunovSolver.prototype.buildEigenstate.Analytic.calcFlux.call(this, matrix, velocity, hTotal, gamma, eTotal);
			EulerEquationGodunovSolver.prototype.buildEigenstate.Analytic.calcEigenvalues.call(this, eigenvalues, velocity, speedOfSound);

			//min eigenvector
			eigenvectors[0][0] = 1;
			eigenvectors[0][1] = velocity - speedOfSound;
			eigenvectors[0][2] = hTotal - speedOfSound * velocity;
			//mid eigenvector
			eigenvectors[1][0] = 1;
			eigenvectors[1][1] = velocity;
			eigenvectors[1][2] = .5 * velocity * velocity;
			//max eigenvector
			eigenvectors[2][0] = 1;
			eigenvectors[2][1] = velocity + speedOfSound;
			eigenvectors[2][2] = hTotal + speedOfSound * velocity;
			
			//calculate eigenvector inverses numerically ... 
			mat33invert(eigenvectorsInverse, eigenvectors);
		}
	},
	Numeric : {
		calcAll : function(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, hTotal, gamma, eTotal) {
			
			//calculate matrix & eigenvalues & vectors at interface from state at interface
			let speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * velocity * velocity));	
			
			//matrix, listed per column
			matrix[0][0] = 0;
			matrix[0][1] = (gamma - 3) / 2 * velocity * velocity;
			matrix[0][2] = velocity * ((gamma - 1) / 2 * velocity * velocity - hTotal);
			matrix[1][0] = 1;
			matrix[1][1] = (3 - gamma) * velocity;
			matrix[1][2] = hTotal - (gamma - 1) * velocity * velocity;
			matrix[2][0] = 0;
			matrix[2][1] = gamma - 1;
			matrix[2][2] = gamma * velocity;

			//analytic eigenvectors and eigenvalues

			//eigenvalues: min, mid, max
			eigenvalues[0] = velocity - speedOfSound;
			eigenvalues[1] = velocity;
			eigenvalues[2] = velocity + speedOfSound;
			//min eigenvector
			eigenvectors[0][0] = 1;
			eigenvectors[0][1] = velocity - speedOfSound;
			eigenvectors[0][2] = hTotal - speedOfSound * velocity;
			//mid eigenvector
			eigenvectors[1][0] = 1;
			eigenvectors[1][1] = velocity;
			eigenvectors[1][2] = .5 * velocity * velocity;
			//max eigenvector
			eigenvectors[2][0] = 1;
			eigenvectors[2][1] = velocity + speedOfSound;
			eigenvectors[2][2] = hTotal + speedOfSound * velocity;		
			//calculate eigenvector inverses ... 
			mat33invert(eigenvectorsInverse, eigenvectors);

			//reconstructed from analytic ... works
			// just don't forget that a_ij = a[j][i] (i.e. transpose your i's and j's and all your other matrix index references)
			/*
			totalError = 0;
			lastReconstructedMatrix = [];
			lastMatrix = [];
			lastError = [];
			for (let i = 0; i < 3; ++i) {
				lastReconstructedMatrix[i] = [];
				lastMatrix[i] = [];
				lastError[i] = [];
				for (let j = 0; j < 3; ++j) {
					let desired = matrix[i][j];
					
					//a_ij = q_ik l_k invq_kj
					let sum = 0;
					for (let k = 0; k < 3; ++k) {
						sum += eigenvectors[k][j] * eigenvalues[k] * eigenvectorsInverse[i][k];
					}
					lastMatrix[i][j] = desired; 
					lastReconstructedMatrix[i][j] = sum;

					lastError[i][j] = Math.abs(sum - desired);
					totalError += lastError[i][j];
				}
			}
			*/
		
			let matrixRowMajor = [];
			for (let i = 0; i < 3; ++i) {
				matrixRowMajor[i] = [];
				for (let j = 0; j < 3; ++j) {
					matrixRowMajor[i][j] = matrix[j][i];
				}
			}
			let u = [];
			let w = [];
			let v = [];

			try {
				//uses row major when my storage is column major
				svd(matrixRowMajor, u/*eigenvectors*/, w/*eigenvalues*/, v/*eigenvectorsInverse*/);
			} catch (e) {
				console.log('failed on matrix',matrix);
				throw e;
			}

/* Uncomment this when you feel safe using this.
 * That may never happen. Turns out HLL solvers don't use eigenvectors! 
 * So numerical decomposition will most likely always be less accurate. */
			for (let i = 0; i < 3; ++i) {
				eigenvalues[i] = w[i];
				for (let j = 0; j < 3; ++j) {
					eigenvectors[i][j] = u[j][i];
					eigenvectorsInverse[i][j] = v[i][j];
				}
			}
			//using v^t causes the system to explode, so 
			//re-calculate eigenvector inverses
			mat33invert(eigenvectorsInverse, eigenvectors);
/**/
			
			/*
			according to this, the difference between my svd-reconstructed matrix
			and the matrix itself is never more than 1e-13 ... 
			so why the oscillations?
			*/
			let totalError = 0;
			let totalEigenError = 0;
			let lastReconstructedMatrix = [];
			let lastMatrix = [];
			let lastError = [];
			let lastEigenvalues = [];
			let lastReconstructedEigenvalues = [];
			for (let i = 0; i < 3; ++i) {
				lastReconstructedMatrix[i] = [];
				lastMatrix[i] = [];
				lastError[i] = [];
				for (let j = 0; j < 3; ++j) {
					let desired = matrix[i][j];
					
					//a_ij = q_ik l_k invq_kj
					let sum = 0;
					for (let k = 0; k < 3; ++k) {
						sum += u[j][k] * w[k] * v[i][k];
					}
					lastMatrix[i][j] = desired; 
					lastReconstructedMatrix[i][j] = sum;

					lastError[i][j] = Math.abs(sum - desired);
					totalError += lastError[i][j];
				}
				lastEigenvalues[i] = eigenvalues[i];
				lastReconstructedEigenvalues[i] = w[i];
				let eigenError = Math.abs(eigenvalues[i] - w[i]);
				totalEigenError += eigenError;
			}
			window.maximalError = Math.max(window.maximalError, totalError);
			window.maximalEigenError = Math.max(window.maximalEigenError, totalEigenError);
			//'desired '+lastMatrix+'\ngot '+lastReconstructedMatrix+'\nerror '+lastError+'\ntotal error '+totalError+' max '+maximalError+'\neig desired '+lastEigenvalues+'\neig got '+lastReconstructedEigenvalues+'\neig total error '+totalEigenError+' max '+maximalEigenError
		}
	}
};
EulerEquationGodunovSolver.prototype.getPrimitives = eulerGetPrimitives;

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
class EulerEquationGodunovExplicit extends EulerEquationGodunovSolver {
	initStep(...args) {
		super.initStep(...args);
		
		for (let ix = 1; ix < this.nx; ++ix) {
			let densityL = this.qR[ix-1][0];
			let velocityL = this.qR[ix-1][1] / densityL;
			let energyTotalL = this.qR[ix-1][2] / densityL;
			let energyKineticL = .5 * velocityL * velocityL;
			let energyInternalL = energyTotalL - energyKineticL;
			let pressureL = (this.gamma - 1) * densityL * energyInternalL;
			let hTotalL = energyTotalL + pressureL / densityL;
			
			let densityR = this.qL[ix][0];
			let velocityR = this.qL[ix][1] / densityR;
			let energyTotalR = this.qL[ix][2] / densityR;
			let energyKineticR = .5 * velocityR * velocityR;
			let energyInternalR = energyTotalR - energyKineticR;
			let pressureR = (this.gamma - 1) * densityR * energyInternalR;
			let hTotalR = energyTotalR + pressureR / densityR;
		
			let velocity = (velocityL + velocityR) * .5;
			let hTotal = (hTotalL + hTotalR) * .5;
			let energyTotal = (energyTotalL + energyTotalR) * .5;
			
			//compute eigenvectors and values at the interface based on averages
			EulerEquationGodunovExplicit.prototype.buildEigenstate[this.eigenDecomposition].calcAll.call(this,
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				velocity, hTotal, this.gamma, energyTotal);
		}	
	}
}

/*
relation:

eigenvalues:
lambda 1 = u - Cs
lambda 2 = u
lambda 3 = u + Cs
eigenvectors:
e 1 = (1, u - Cs, h_total - Cs u)
e 2 = (1, u, .5*u^2)
e 3 = (1, u + Cs, h_total + Cs u)
rho = q0
u = q1 / q0
h_total = e_total + P / rho
Cs = sqrt(gamma P / rho)
*/
class EulerEquationRoeExplicit extends EulerEquationGodunovSolver {
	initStep(...args) {
		super.initStep(...args);
		for (let ix = 1; ix < this.nx; ++ix) {
			//compute Roe averaged interface values
			const densityL = this.qR[ix-1][0];
			const velocityL = this.qR[ix-1][1] / densityL;
			const energyTotalL = this.qR[ix-1][2] / densityL;
			const energyKineticL = .5 * velocityL * velocityL;
			const energyInternalL = energyTotalL - energyKineticL;
			const pressureL = (this.gamma - 1) * densityL * energyInternalL;
			const hTotalL = energyTotalL + pressureL / densityL;
			const roeWeightL = Math.sqrt(densityL);
			
			const densityR = this.qL[ix][0];
			const velocityR = this.qL[ix][1] / densityR;
			const energyTotalR = this.qL[ix][2] / densityR;
			const energyKineticR = .5 * velocityR * velocityR;
			const energyInternalR = energyTotalR - energyKineticR;
			const pressureR = (this.gamma - 1) * densityR * energyInternalR;
			const hTotalR = energyTotalR + pressureR / densityR;
			const roeWeightR = Math.sqrt(densityR);
			
			const denom = roeWeightL + roeWeightR;
			const invDenom = 1 / denom;
			const velocity = (roeWeightL * velocityL + roeWeightR * velocityR) * invDenom;
			const hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) * invDenom;
			const energyTotal = (roeWeightL * energyTotalL + roeWeightR * energyTotalR) * invDenom;

			//compute eigenvectors and values at the interface based on Roe averages
			EulerEquationRoeExplicit.prototype.buildEigenstate[this.eigenDecomposition].calcAll.call(this,
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				velocity, hTotal, this.gamma, energyTotal);
		}
	}
}

class EulerEquationHLLExplicit extends EulerEquationGodunovSolver {
	initStep(...args) {
		super.initStep(...args);
		
		for (let ix = 2; ix < this.nx - 1; ++ix) {
			//compute Roe averaged interface values
			const densityL = this.qR[ix-1][0];
			const velocityL = this.qR[ix-1][1] / densityL;
			const energyTotalL = this.qR[ix-1][2] / densityL;
			const energyKineticL = .5 * velocityL * velocityL;
			const energyInternalL = energyTotalL - energyKineticL;
			const pressureL = (this.gamma - 1) * densityL * energyInternalL;
			const hTotalL = energyTotalL + pressureL / densityL;
			const roeWeightL = Math.sqrt(densityL);
			
			const densityR = this.qL[ix][0];
			const velocityR = this.qL[ix][1] / densityR;
			const energyTotalR = this.qL[ix][2] / densityR;
			const energyKineticR = .5 * velocityR * velocityR;
			const energyInternalR = energyTotalR - energyKineticR;
			const pressureR = (this.gamma - 1) * densityR * energyInternalR;
			const hTotalR = energyTotalR + pressureR / densityR;
			const roeWeightR = Math.sqrt(densityR);
			
			const denom = roeWeightL + roeWeightR;
			const invDenom = 1 / denom;
			const velocity = (roeWeightL * velocityL + roeWeightR * velocityR) * invDenom;
			const hTotal = (roeWeightL * hTotalL + roeWeightR * hTotalR) * invDenom;
			const speedOfSound = Math.sqrt((this.gamma - 1) * (hTotal - .5 * velocity * velocity));

			//used by superclass for determining timestep
			//( should I be using the Roe wavespeeds or should I be using the HLL flux ... divided by the state variable (Hydrodynamics II 7.44) to determine CFL?
			this.interfaceEigenvalues[ix][0] = velocity - speedOfSound;
			this.interfaceEigenvalues[ix][1] = velocity;
			this.interfaceEigenvalues[ix][2] = velocity + speedOfSound;
		}
	}
	step(dt) {
		const deriv = EulerEquationHLLExplicit.prototype.calcDerivative;
		explicitMethods[this.explicitMethod].call(this, dt, deriv);
	}
	calcDerivative(dt, dq_dt) {
		const nx = this.nx;
		for (let i = 2; i < nx - 1; ++i) {
			const qL = this.q[i-1];
			const densityL = qL[0];
			const velocityL = qL[1] / densityL;
			const velocitySqL = velocityL * velocityL;
			const energyTotalL = qL[2] / densityL;
			const energyKineticL = .5 * velocitySqL;
			const energyInternalL = energyTotalL - energyKineticL;
			const pressureL = (this.gamma - 1) * densityL * energyInternalL;
			const hTotalL = energyTotalL + pressureL / densityL;
			const roeWeightL = Math.sqrt(densityL);
			const speedOfSoundL = Math.sqrt((this.gamma - 1) * (hTotalL - energyKineticL));
			const eigenvalueLMin = velocityL - speedOfSoundL;
			const eigenvalueLMax = velocityL + speedOfSoundL;

			const fluxL = [
				densityL * velocityL,
				densityL * velocityL * velocityL + pressureL,
				densityL * velocityL * hTotalL
			];

			const qR = this.q[i];
			const densityR = qR[0];
			const velocityR = qR[1] / densityR;
			const velocitySqR = velocityR * velocityR;
			const energyTotalR = qR[2] / densityR;
			const energyKineticR = .5 * velocitySqR;
			const energyInternalR = energyTotalR - energyKineticR;
			const pressureR = (this.gamma - 1) * densityR * energyInternalR;
			const hTotalR = energyTotalR + pressureR / densityR;
			const roeWeightR = Math.sqrt(densityR);
			const speedOfSoundR = Math.sqrt((this.gamma - 1) * (hTotalR - energyKineticR));
			const eigenvalueRMin = velocityR - speedOfSoundR;
			const eigenvalueRMax = velocityR + speedOfSoundR;

			const fluxR = [
				densityR * velocityR,
				densityR * velocityR * velocityR + pressureR,
				densityR * velocityR * hTotalR
			];

			const invDenom = 1 / (roeWeightL + roeWeightR);
			const velocityC = (roeWeightL * velocityL + roeWeightR * velocityR) * invDenom;
			const hTotalC = (roeWeightL * hTotalL + roeWeightR * hTotalR) * invDenom;
			const speedOfSoundC = Math.sqrt((this.gamma - 1) * (hTotalC - .5 * velocityC * velocityC));

			const eigenvalueCMin = velocityC - speedOfSoundC;
			const eigenvalueCMax = velocityC + speedOfSoundC;

			const sL = Math.min(eigenvalueLMin, eigenvalueCMin);
			const sR = Math.max(eigenvalueRMax, eigenvalueCMax);

			for (let j = 0; j < 3; ++j) {
				if (sL >= 0) {
					this.interfaceFlux[i][j] = fluxL[j];
				} else if (sL <= 0 && sR >= 0) {
					this.interfaceFlux[i][j] = (sR * fluxL[j] - sL * fluxR[j] + sL * sR * (qR[j] - qL[j])) / (sR - sL);
				} else if (sR <= 0) {
					this.interfaceFlux[i][j] = fluxR[j];
				}
			}

			//TODO flux limiter
			// so how do we deconstruct the flux change across the interface such that we can compute its ratios and apply our limiter ...
		}

		//zero boundary flux
		for (let j = 0; j < 3; ++j) {
			this.interfaceFlux[0][j] = this.interfaceFlux[this.nx][j] = 0;
		}

		//update cells
		for (let i = this.nghost; i < this.nx-this.nghost; ++i) {
			for (let j = 0; j < 3; ++j) {
				dq_dt[i][j] = -(this.interfaceFlux[i+1][j] - this.interfaceFlux[i][j]) / (this.xi[i+1] - this.xi[i]);
			}
		}
		for (let i = 0; i < this.nghost; ++i) {
			for (let j = 0; j < 3; ++j) {
				dq_dt[i][j] = 0;
				dq_dt[this.nx-i-1][j] = 0;
			}
		}
	}
}

let hdSimulation = {
	methods : {
		'Burgers / Explicit' : EulerEquationBurgersExplicit.prototype,
		'Burgers / Backward Euler via Gauss Seidel' : EulerEquationBurgersBackwardEulerGaussSeidel.prototype,
		//'Burgers / Backward Euler via Block Tridiagonal' : EulerEquationBurgersBackwardEulerTridiagonal.prototype,
		//'Godunov / Explicit' : EulerEquationGodunovExplicit.prototype,
		'HLL / Explicit' : EulerEquationHLLExplicit.prototype,
		'Roe / Explicit' : EulerEquationRoeExplicit.prototype
	},
	initialConditions : {
		Sod : function() {
			this.resetCoordinates(-1, 1);
			for (let i = 0; i < this.nx; ++i) {
				const x = this.x[i];
				const rho = (x < (xmin * .7 + xmax * .3)) ? 1 : .1;
				const u = 0;
				const eTotal = 1;
				this.q[i][0] = rho; 
				this.q[i][1] = rho * u; 
				this.q[i][2] = rho * eTotal; 
			}
		},
		Sedov : function() {
			this.resetCoordinates(-1, 1);
			const pressure = 1e-5;
			const rho = 1;
			for (let i = 0; i < this.nx; ++i) {
				this.q[i][0] = rho;
				this.q[i][1] = 0;
				this.q[i][2] = pressure / (this.gamma - 1);
			}
			this.q[Math.floor(this.nx/2)][2] = 1e+5;
		},
		Advect : function() {
			this.resetCoordinates(-1, 1);
			const xmid = .5 * (xmax + xmin);
			for (let i = 0; i < this.nx; ++i) {
				const x = this.x[i];
				const xGreaterThanMid = x > xmid;
				const rho = xGreaterThanMid ? .5 : 1;
				const u = 1;
				const energyKinetic = .5 * u * u;
				const pressure = 1;
				const energyTotal = pressure / (rho * (this.gamma - 1)) + energyKinetic;
				this.q[i][0] = rho;
				this.q[i][1] = rho * u;
				this.q[i][2] = rho * energyTotal;
			}
		},
		Wave : function() {
			this.resetCoordinates(-1, 1);
			const xmid = .5 * (xmin + xmax);
			const dg = .1 * (xmax - xmin);
			for (let i = 0; i < this.nx; ++i) {
				const x = this.x[i];
				const dx = x - xmid;
				const rho = 1 + .3 * Math.exp(-(dx*dx)/(dg*dg));
				const u = 0;
				const eTotal = 1;
				this.q[i][0] = rho;
				this.q[i][1] = rho * u;
				this.q[i][2] = rho * eTotal;
			}
		}
	}
};

/*
From "Numerical Relativity", Alcubierre. 2008. ch.10 "Examples of Numerical Relativity" section 2: "Toy 1+1 Relativity"

State variables:
q[0] = D_alpha = d/dx ln alpha
q[1] = D_g = d/dx ln g
q[2] = K^tilde = sqrt(g) K

... for primtiive variables
alpha = lapse
g = g_xx = metric
K = extrinsic curvature = K^x_x

...and some extra variables:
f comes from Bona-Masso slicing family: d/dt alpha = -alpha^2 f(alpha) K
(shift value is zero)
*/

const adm_BonaMasso_f = 1; 		//aiming to recreate the f=.5, f=1, f=1.5 graphs ...

function admEquationsBuildEigenstate(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, alpha, g) {
	const f = adm_BonaMasso_f;
	
	const oneOverSqrtG = 1 / Math.sqrt(g);
	const sqrtF = Math.sqrt(f);
	
	matrix[0][0] = 0;
	matrix[0][1] = 0;
	matrix[0][2] = alpha * oneOverSqrtG;
	matrix[1][0] = 0;
	matrix[1][1] = 0;
	matrix[1][2] = 0;
	matrix[2][0] = alpha * f * oneOverSqrtG;
	matrix[2][1] = 2 * alpha * oneOverSqrtG;
	matrix[2][2] = 0;

	eigenvalues[0] = -alpha * sqrtF * oneOverSqrtG;
	eigenvalues[1] = 0;
	eigenvalues[2] = -eigenvalues[0];

	//column 0: min
	eigenvectors[0][0] = f;
	eigenvectors[0][1] = 2;
	eigenvectors[0][2] = -sqrtF;

	//column 1: mid
	eigenvectors[1][0] = 0;
	eigenvectors[1][1] = 1;
	eigenvectors[1][2] = 0;

	//column 2: max
	eigenvectors[2][0] = f;
	eigenvectors[2][1] = 2;
	eigenvectors[2][2] = sqrtF;
	
	//calculate eigenvector inverses ... 
	mat33invert(eigenvectorsInverse, eigenvectors);
}

class ADMRoeExplicit extends GodunovSolver {
	initStep(...args) {
		super.initStep(...args);
		const dx = (xmax - xmin) / this.nx;
		for (let ix = 2; ix < this.nx-1; ++ix) {
			const ixL = ix-1;
			const ixR = ix;

			//q_ix,0 = d/dx ln alpha
			const ln_alpha_L = (this.q[ixL+1][0] - this.q[ixL-1][0]) / (2 * dx);
			const alpha_L = Math.exp(ln_alpha_L);

			//q_ix,1 = d/dx ln g
			const ln_g_L = (this.q[ixL+1][1] - this.q[ixL-1][1]) / (2 * dx);
			const g_L = Math.exp(ln_g_L);

			const ln_alpha_R = (this.q[ixR+1][0] - this.q[ixR-1][0]) / (2 * dx);
			const alpha_R = Math.exp(ln_alpha_R);
			
			const ln_g_R = (this.q[ixR+1][1] - this.q[ixR-1][1]) / (2 * dx);
			const g_R = Math.exp(ln_g_R);

			//TODO pick a better weighting value
			const alpha = .5 * (alpha_L + alpha_R);
			const g = .5 * (g_L + g_R);

			//compute eigenvectors and values at the interface based on Roe averages
			admEquationsBuildEigenstate(
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				alpha, g);
		}
		//how about those boundary eigenstates?
	}
	getPrimitives(i) {
		const dx = (xmax - xmin) / this.nx;
		const nx = this.nx;
		
		const iL = i <= 0 ? 0 : i - 1;
		const iR = i >= nx-1 ? nx-1 : i + 1;
		
		const ln_alpha = (this.q[iR][0] - this.q[iL][0]) / (2 * dx);
		const alpha = Math.exp(ln_alpha);

		//q_i,1 = d/dx ln g
		const ln_g = (this.q[iR][1] - this.q[iL][1]) / (2 * dx);
		const g = Math.exp(ln_g);

		const KTilde = this.q[i][2];
		const K = KTilde / Math.sqrt(g);
		
		return [alpha, g, K];
	}
}

const admSimulation = {
	methods : {
		'Roe / Explicit' : ADMRoeExplicit.prototype
	},
	initialConditions : {
		GaugeShock : function() {
			this.resetCoordinates(-30, 30);
			const xmid = (xmax + xmin) * .5;
			for (let i = 0; i < this.nx; ++i) {
				const x = (this.x[i] - xmid) / ((xmax - xmid) / 3);
				const h = Math.exp(-x*x); 
				const dh_dx = -2 * x * h;
				const d2h_dx2 = 2 * h * (2 * x * x - 1);
				const g = 1 - dh_dx * dh_dx;
				const D_g = -2 * dh_dx * d2h_dx2 / g;
				const KTilde = -d2h_dx2 / g;
				const f = adm_BonaMasso_f;
				const D_alpha = Math.sqrt(f) * KTilde;
				this.q[i][0] = D_alpha;
				this.q[i][1] = D_g;
				this.q[i][2] = KTilde;
			}
		}
	}
};

const mu0 = 1;	//permittivity of free space
const sqrtMu0 = Math.sqrt(mu0);

function mhdBuildEigenstate(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, alpha, g) {
	
	matrix[0][0] = 0;
	matrix[0][1] = 0;
	matrix[0][2] = alpha * oneOverSqrtG;
	
	matrix[1][0] = 0;
	matrix[1][1] = 0;
	matrix[1][2] = 0;
	
	matrix[2][0] = alpha * f * oneOverSqrtG;
	matrix[2][1] = 2 * alpha * oneOverSqrtG;
	matrix[2][2] = 0;

	eigenvalues[0] = -alpha * sqrtF * oneOverSqrtG;
	eigenvalues[1] = 0;
	eigenvalues[2] = -eigenvalues[0];

	//column 0: min
	eigenvectors[0][0] = f;
	eigenvectors[0][1] = 2;
	eigenvectors[0][2] = -sqrtF;

	//column 1: mid
	eigenvectors[1][0] = 0;
	eigenvectors[1][1] = 1;
	eigenvectors[1][2] = 0;

	//column 2: max
	eigenvectors[2][0] = f;
	eigenvectors[2][1] = 2;
	eigenvectors[2][2] = sqrtF;

	//...and the inverse...
}


class MHDGodunovExplicit extends GodunovSolver {
	initStep(...args) {
		super.initStep(...args);
		const dx = (xmax - xmin) / this.nx;
		//same idea as Godunov but with Roe weighting: sqrt(rho)
		for (let ix = 2; ix < this.nx-1; ++ix) {
			const ixL = ix-1;
			const ixR = ix;

			//q_ix,0 = d/dx ln alpha
			const ln_alpha_L = (this.q[ixL+1][0] - this.q[ixL-1][0]) / (2 * dx);
			const alpha_L = Math.exp(ln_alpha_L);

			//q_ix,1 = d/dx ln g
			const ln_g_L = (this.q[ixL+1][0] - this.q[ixL-1][0]) / (2 * dx);
			const g_L = Math.exp(ln_g_L);

			const ln_alpha_R = (this.q[ixR+1][0] - this.q[ixR-1][0]) / (2 * dx);
			const alpha_R = Math.exp(ln_alpha_R);
			
			const ln_g_R = (this.q[ixR+1][0] - this.q[ixR-1][0]) / (2 * dx);
			const g_R = Math.exp(ln_g_R);

			const alpha = .5 * (alpha_L + alpha_R);
			const g = .5 * (g_L + g_R);

			//compute eigenvectors and values at the interface based on Roe averages
			mhdBuildEigenstate(
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix], 
				alpha, g);
		}
		//how about those boundary eigenstates?
	}
	step(dt) {
		const deriv = MHDGodunovExplicit.prototype.calcDerivative;
		explicitMethods[this.explicitMethod].call(this, dt, deriv);
	}
	getPrimitives(i) {
		const rho = this.q[0];
		return [
			rho, this.q[1] / rho, this.q[2] / rho, this.q[3] / rho,
			this.q[4] / rho, this.q[5] * sqrtMu0, this.q[6] * sqrtMu0, this.q[7] * sqrtMu0
		];
	}
}

let mhdSimulation = {
	methods : {
		'Godunov / Explicit' : MHDGodunovExplicit.prototype
	},
	initialConditions : {
		Sod : function() {
			this.resetCoordinates(-1, 1);
			for (let i = 0; i < this.nx; ++i) {
				const x = this.x[i];
				const rho = (x < (xmin * .7 + xmax * .3)) ? 1 : .1;
				const u = 0;
				const v = 0;
				const w = 0;
				const eTotal = 1;
				const Bx = 0;
				const By = 0;
				const Bz = 0;
				this.q[i][0] = rho; 
				this.q[i][1] = rho * u; 
				this.q[i][2] = rho * v;
				this.q[i][3] = rho * w;
				this.q[i][4] = rho * eTotal; 
				this.q[i][5] = Bx / sqrtMu0;
				this.q[i][6] = By / sqrtMu0;
				this.q[i][7] = Bz / sqrtMu0;
			}
		}
	}
};

//"A Flux-Split Algorithm applied to Relativistic Flows",  R. Donat, J. A. Font, J. M. Ibanez, A. Marquina
let srhdBuildEigenstate = function(matrix, eigenvalues, eigenvectors, eigenvectorsInverse, velocity, W, density, energyInternal, pressure, gamma) {

	//for newtonian, for enthalpy based on energy total:
	//let hTotal = energyTotal + pressure / density;
	//let speedOfSound = Math.sqrt((gamma - 1) * (hTotal - .5 * velocity * velocity));	

	const enthalpy = 1 + energyInternal + pressure / density;
	//Alcubierre 7.6.53
	const speedOfSound = Math.sqrt(gamma * (gamma - 1) * energyInternal / (1 + gamma * energyInternal));

	const denom = 1 - velocity * velocity * speedOfSound * speedOfSound;
	const a = velocity * (1 - speedOfSound * speedOfSound);
	const b = speedOfSound * (1 - velocity * velocity);

	//eigenvalues: min, mid, max
	eigenvalues[0] = (a - b) / denom;
	eigenvalues[1] = velocity;
	eigenvalues[2] = (a + b) / denom; 

	//const kappaTilde = kappa / density;
	const K = 1;//kappaTilde / (kappaTilde - speedOfSound * speedOfSound);
	const Aminus= (1 - velocity * velocity) / (1 - velocity * eigenvalues[0]);
	const Aplus = (1 - velocity * velocity) / (1 - velocity * eigenvalues[2]);
	
	//min eigenvector
	eigenvectors[0][0] = 1;
	eigenvectors[0][1] = enthalpy * W * Aminus * eigenvalues[0]; 
	eigenvectors[0][2] = enthalpy * W * Aminus - 1;
	//mid eigenvector
	eigenvectors[1][0] = K / (enthalpy * W);
	eigenvectors[1][1] = velocity;
	eigenvectors[1][2] = 1 - K / (enthalpy * W);
	//max eigenvector
	eigenvectors[2][0] = 1;
	eigenvectors[2][1] = enthalpy * W * Aplus * eigenvalues[2]; 
	eigenvectors[2][2] = enthalpy * W * Aplus - 1; 

/* manual -- producing NaNs * /
	//calculate eigenvector inverses numerically ... though I don't have to ...
	//mat33invert(eigenvectorsInverse, eigenvectors);
	//(singular) or we can do it symbolically...
	const delta = enthalpy * enthalpy * enthalpy * W * (K - 1) * (1 - velocity * velocity) * (Aplus * eigenvalues[2] - Aminus * eigenvalues[0]);
	const ls = enthalpy * enthalpy / delta;
	//min
	eigenvectorsInverse[0][0] = -ls * (enthalpy * W * Aminus * (velocity - eigenvalues[0]) - velocity + k * Aminus * eigenvalues[0]);
	eigenvectorsInverse[0][1] = -ls * (1 - K * Aminus);
	eigenvectorsInverse[0][2] = -ls * (-velocity + K * Aminus * eigenvalues[0]);
	//mid
	const lms = W / (K - 1);
	eigenvectorsInverse[1][0] = lms * (enthalpy - W);
	eigenvectorsInverse[1][1] = lms * W * velocity;
	eigenvectorsInverse[1][2] = lms * -W;
	//max
	eigenvectorsInverse[2][0] = ls * (enthalpy * W * Aplus * (velocity - eigenvalues[2]) - velocity + k * Aplus * eigenvalues[2]);
	eigenvectorsInverse[2][1] = ls * (1 - K * Aplus);
	eigenvectorsInverse[2][2] = ls * (-velocity + K * Aplus * eigenvalues[2]);
/**/
/* numeric - producing singular values */
	mat33invert(eigenvectorsInverse, eigenvectors);
/**/
	//calculate flux numerically ... though I don't have to ...
	for (let i = 0; i < 3; ++i) {
		for (let j = 0; j < 3; ++j) {
			//a_ij = u_ik w_k vt_kj
			let s = 0;
			for (let k = 0; k < 3; ++k) {
				s += eigenvectors[k][i] * eigenvalues[k] * eigenvectorsInverse[j][k];
			}
			matrix[j][i] = s;
		}
	}
}

class SRHDGodunovExplicit extends GodunovSolver {
	initStep(...args) {
		super.initStep(...args);
		const dx = (xmax - xmin) / this.nx;
		
		const prims = [];
		for (let i = 1; i < this.nx-1; ++i) {
			prims[i] = this.getPrimitives(i);
		}
		
		for (let ix = 2; ix < this.nx-1; ++ix) {
			const ixL = ix-1;
			const ixR = ix;

			const primsL = prims[ixL];
			const primsR = prims[ixR];
			const primsAvg = {};
			for (let k in primsL) {
				primsAvg[k] = (primsL[k] + primsR[k]) / 2;
			}

			srhdBuildEigenstate(
				this.interfaceMatrix[ix],
				this.interfaceEigenvalues[ix], 
				this.interfaceEigenvectors[ix], 
				this.interfaceEigenvectorsInverse[ix],
				primsAvg.v, primsAvg.W, primsAvg.rho0, primsAvg.eInternal, primsAvg.p, this.gamma);
		}
		//how about those boundary eigenstates?
	}
	step(dt) {
		let deriv = SRHDGodunovExplicit.prototype.calcDerivative;
		explicitMethods[this.explicitMethod].call(this, dt, deriv);
	}
	getPrimitives(i) {
		//actually a nonlinear problem ... so ... solve this one later
		const q = this.q[i];
		const D = q[0];
		const S = q[1];
		const E = q[2];
		const p = q.pressure;	//last iteration
		const v = S / (E + D + p);
		const W = 1 / Math.sqrt(1 - v * v);
		const rho0 = D / W;
		const energyInternal = (E + D * (1 - W) + p * (1 - W * W)) / (D * W);
		const nextP = (this.gamma - 1) * rho0 * energyInternal;
		//rinse, lather, repeat
		p = nextP;
		q.pressure = p;
		const prims = {};
		prims.v = v;
		prims.W = W;
		prims.rho0 = rho0;
		prims.eInternal = energyInternal;
		prims.p = p;
		//for graphing...
		prims[0] = rho0;
		prims[1] = v;
		prims[2] = energyInternal;
		return prims;
	}
}


/*
Using Alcubierre, Baumgarte & Shapiro, and...
http://relativity.livingreviews.org/Articles/lrr-2003-7/download/lrr-2003-7Color.pdf
*/
let srhdSimulation = {
	methods : {
		'Godunov / Explicit' : SRHDGodunovExplicit.prototype
	},
	initialConditions : {
		Sod : function() {
			this.resetCoordinates(-1, 1);
			for (let i = 0; i < this.nx; ++i) {
				const x = this.x[i];
				const rho0 = (x < (xmin * .7 + xmax * .3)) ? 1 : .1;
				const v = 0;	//three-velocity, in geometric units, so vx = 1 <=> speed of light
				const W = 1 / Math.sqrt(1 - v * v);
			
				const eInternal = 1;
				const p = (this.gamma - 1) * rho0 * eInternal;	//
				const h = 1 + eInternal + p / rho0;
				
				const D = rho0 * W; 
				const S = rho0 * h * W * W * v; 
				const E = rho0 * h * W * W - p - D;;
				this.q[i][0] = D;
				this.q[i][1] = S;
				this.q[i][2] = E;
				this.q[i].pressure = p;
			}
		}
	}
};

//hmm, maybe I should combine all of these into one list, and have them individuall state who they belong to
let simulations = {
	Euler : hdSimulation
	//MHD : mhdSimulation,
	//SRHD : srhdSimulation,
	//ADM : admSimulation
};

class HydroState {
	//call this to reflect this.x and this.xi after changing xmin, xmax, or this.nx
	resetCoordinates(xmin_, xmax_) {
		xmin = xmin_;
		xmax = xmax_;
		onresize();

		//x_i: cell positions
		this.x = new Float32Array(this.nx);
		for (let i = 0; i < this.nx; ++i) {
			this.x[i] = xmin + (xmax - xmin) * i / (this.nx-1);
		}
		
		//x_{i-1/2}: interface positions
		this.xi = new Float32Array(this.nx+1);
		for (let i = 1; i < this.nx+1; ++i) {
			this.xi[i] = .5*(this.x[i] + this.x[i-1]);
		}
		this.xi[0] = 2 * this.xi[1] - this.xi[2];
		this.xi[this.nx] = 2 * this.xi[this.nx-1] - this.xi[this.nx-2]; 
	}
	constructor(args) {
		this.nx = args.size;
		this.cfl =.5;
		
		//used for Euler equation of state
		this.gamma = args.gamma;
	
		this.resetCoordinates();

		//q_j,i: state vector, stored as q[j][i]
		//q_0,i: density: rho
		//q_1,i: momentum: rho * v
		//q_2,i: work: rho * e
		this.q = [];
		for (let i = 0; i < this.nx; ++i) {
			this.q[i] = [];
		}

		this.tmpq = [];
		for (let j = 0; j < 5; ++j) {
			this.tmpq[j] = [];
			for (let i = 0; i < this.nx; ++i) {
				this.tmpq[j][i] = [];
			}
		}


		hdSimulation.initialConditions.Sod.call(this);
		
		//p_i: pressure
		this.pressure = new Float32Array(this.nx);
	
		
		//used for Burgers
		
		
		//r_{i-1/2}	
		this.r = [];
		for (let i = 0; i < this.nx+1; ++i) {
			this.r[i] = [0,0,0];
		}
		
		//f_{i-1/2}: cell flux
		this.interfaceFlux = [];
		for (let i = 0; i < this.nx+1; ++i) {
			this.interfaceFlux[i] = [0,0,0];
		}
	
		//used with Burgers 
		//u_{i-1/2}: interface velocity
		this.ui = new Float32Array(this.nx+1);
	
		//used with Godunov
		this.qMid = [];
		for (let i = 0; i <= this.nx; ++i) this.qMid[i] = [0,0,0];
		//q_i,L and q_i,R	<- L and R at cell (as PPM uses), not at interface (as Godunov uses)
		//	so without interpolation, q_{i-1/2},R = q_i,L and q_{i-1/2},L = q_{i-1},R
		this.qL = [];
		for (let i = 0; i < this.nx; ++i) this.qL[i] = [0,0,0];
		this.qR = [];
		for (let i = 0; i < this.nx; ++i) this.qR[i] = [0,0,0];

		//only used with Godunov and Roe eqn advection code:
		//calculated before dt
		// used for dt and for fluxTilde
		this.interfaceMatrix = [];
		this.interfaceEigenvalues = [];	//lambda_{i-1/2},j: interface eigenvalues
		this.interfaceEigenvectors = [];		//e_{i-1/2},j,k: interface eigenvectors
		this.interfaceEigenvectorsInverse = [];		//[e_{i-1/2},j,k]^-1: interface eigenvector column matrix inverse 
		for (let ix = 0; ix < this.nx+1; ++ix) {
			this.interfaceMatrix[ix] = [[1,0,0], [0,1,0], [0,0,1]];
			this.interfaceEigenvalues[ix] = [0, 0, 0];
			this.interfaceEigenvectors[ix] = [[1,0,0], [0,1,0], [0,0,1]];
			this.interfaceEigenvectorsInverse[ix] = [[1,0,0], [0,1,0], [0,0,1]];
		}

		//used for Godunov
		this.rTilde = [];
		for (let i = 0; i < this.nx+1; ++i) {
			this.rTilde[i] = [0,0,0];
		}

		//state change over interface in interface Roe eigenspace
		//used for Godunov
		//tilde means in basis of eigenvectors
		this.interfaceDeltaQTilde = [];
		for (let i = 0; i <= this.nx; ++i) {
			this.interfaceDeltaQTilde[i] = [0,0,0];
		}

	
		//used for Backward Euler + Gauss Seidel
		this.oldQ = [];
		for (let i = 0; i < this.nx; ++i) {
			this.oldQ[i] = [0, 0, 0];
		}
		
		//number of ghost cells
		this.nghost = 2;

		//solver configuration
		this.boundaryMethod = 'mirror';
		this.fluxMethod = 'superbee';
		this.simulation = 'Euler';
		this.algorithm = 'Roe / Explicit';
		this.eigenDecomposition = 'Analytic';
		this.explicitMethod = 'RK4';
	}
	boundary() {
		boundaryMethods[this.boundaryMethod](this.nx, this.q);
	}
	step(dt) {
		//apply boundary conditions
		this.boundary();
	
		//solve
		simulations[this.simulation].methods[this.algorithm].step.call(this, dt);
	}
	getPrimitives(...args) {
		return simulations[this.simulation].methods[this.algorithm].getPrimitives.apply(this, args);
	}
	update() {
		//do any pre-calcCFLTimestep preparation (Roe computes eigenvalues here)
		simulations[this.simulation].methods[this.algorithm].initStep.call(this);
		
		//get timestep
		let dt;
		if (useCFL) {
			dt = simulations[this.simulation].methods[this.algorithm].calcCFLTimestep.call(this);
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
		if (!isFinite(size) || !size) size = 200;
		let gamma = Number(urlparams.get('gamma'));
		if (!gamma || !isFinite(gamma)) gamma = 7/5;

		this.state = new HydroState({
			size : size,
			gamma : gamma
		});
		
		//geometry buffers
		this.vertexXBuffer = new glutil.ArrayBuffer({
			dim : 1,
			count : size,
			usage : gl.DYNAMIC_DRAW
		});

		let colors = [[1,0,0,1], [0,1,0,1], [0,.5,1,1], [1,1,1,1]];

		this.primitiveBuffers = [];
		this.stateGraphObjs = [];
		for (let i = 0; i < 3; ++i) {
			let stateBuffer = new glutil.ArrayBuffer({
				dim : 1,
				count : size,
				usage : gl.DYNAMIC_DRAW
			});
			this.primitiveBuffers.push(stateBuffer);
	
			//make graphs
			let stateGraphObj = new glutil.SceneObject({
				mode : gl.LINE_STRIP,
				attrs : {
					vertex : this.vertexXBuffer,
					state : stateBuffer
				},
				shader : graphShader,
				uniforms : {
					color : colors[i % colors.length]
				},
				blend : [gl.ONE, gl.ONE]
			});
			this.stateGraphObjs.push(stateGraphObj);
		}
	}
	update() {
		//update a copy of the grid and its once-refined
		//...and a once-unrefined ... over mergeable cells only?
		//then test for errors and split when needed
		this.state.update();
	
		//update geometry
		for (let i = 0; i < this.state.nx; ++i) {
			this.vertexXBuffer.data[i] = this.state.x[i];
			//rescale to -1,1
			this.vertexXBuffer.data[i] -= xmin;
			this.vertexXBuffer.data[i] *= 2 / (xmax - xmin);
			this.vertexXBuffer.data[i]--;
		}
		this.vertexXBuffer.updateData();
	
		for (let i = 0; i < this.stateGraphObjs.length; ++i) {
			this.stateGraphObjs[i].uniforms.scale = 1/(xmax - xmin);
		}

		//TODO this is only for Euler -- make it generic
		for (let i = 0; i < this.state.nx; ++i) {
			let prims = this.state.getPrimitives(i);
			this.primitiveBuffers[0].data[i] = prims[0];
			this.primitiveBuffers[1].data[i] = prims[1];
			this.primitiveBuffers[2].data[i] = prims[2];
		}
		for (let j = 0; j < 3; ++j) {
			this.primitiveBuffers[j].updateData();
		}
	}
}

let hydro;

function update() {
	//iterate
	if (!pause) hydro.update();
	//draw
	glutil.draw();
	requestAnimationFrame(update);
}

function onresize() {
	glutil.canvas.width = window.innerWidth;
	glutil.canvas.height = window.innerHeight;
	ids.content.style.height = (window.innerHeight - 50)+'px';
	glutil.resize();
}

function buildSelect(id, key, map) {
	let select = ids[id];
	for (let k in map) {
		let option = Option({innerText : k, appendTo : select});
		if (hydro.state[key] == k) {
			option.setAttribute('selected', 'true');
		}
	}
	select.addEventListener('change', e => {
		hydro.state[key] = select.value;
	});
}

canvas = Canvas({
	style : {
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
	removeFromParent(ids.panel);
	removeFromParent(canvas);
	show(ids.webglfail);
	throw e;
}

glutil.view.ortho = true;
glutil.view.zNear = -1;
glutil.view.zFar = 1;
glutil.view.pos[0] = 0;//(xmax + xmin) / 2;
glutil.view.pos[1] = 0;//(ymax + ymin) / 2;
glutil.view.fovY = ymax - ymin;

plainShader = new glutil.Program({
	vertexCode : `
in vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex.xy, 0., 1.);
	gl_PointSize = 3.;
}
`,
	fragmentCode : `
uniform vec4 color;
out vec4 fragColor;
void main() {
	fragColor = color;
}
`,
	uniforms : {
		color : [1,1,1,1]
	}
});

graphShader = new glutil.Program({
	vertexCode : `
in float vertex;
in float state;
uniform float scale;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	gl_Position = projMat * mvMat * vec4(vertex, scale * state, 0., 1.);
	gl_PointSize = 3.;
}
`,
	fragmentCode : `
uniform vec4 color;
out vec4 fragColor;
void main() {
	fragColor = color;
}
`,
	uniforms : {
		color : [1,1,1,1],
		scale : 1
	}
});

hydro = new Hydro();
window.hydro = hydro;

panel = ids.panel;	

let panelContent = ids.content;
ids.menu.addEventListener('click', e => {
	if (hidden(panelContent)) {
		show(panelContent);
	} else {
		hide(panelContent);
	}
});

ids.pause.addEventListener('click', e => {
	pause = !pause;
	ids.pause.setAttribute('src', pause ? 'play.png' : 'pause.png');
});

for (let simulationName in simulations) {
	let parent = ids[simulationName+'_controls'];
	let simulation = simulations[simulationName];
	for (let initialConditionName in simulation.initialConditions) {
		let method = simulation.initialConditions[initialConditionName];
		Button({
			innerText : 'Reset '+initialConditionName,
			events : {
				click : e => {
					method.call(hydro.state);
				},
			},
			appendTo : parent,
		});
		Br({appendTo : parent});
	}
}

let select = ids.gridsize;
[20, 50, 100, 200, 500, 1000].forEach(gridsize => {
	let option = Option({innerText : gridsize, appendTo : select});
	if (hydro.state.nx == gridsize) {
		option.setAttribute('selected', 'true');
	}
});
select.addEventListener('change', e => {
	const params = new URLSearchParams(urlparams);
	params.set('size', select.value);
	location.href = location.origin + location.pathname + '?' + params.toString();
});

buildSelect('boundaryMethod', 'boundaryMethod', boundaryMethods);
buildSelect('fluxMethod', 'fluxMethod', fluxMethods);

buildSelect('Euler_algorithm', 'algorithm', hdSimulation.methods);	//this will change as 'simulations' changes
buildSelect('SRHD_algorithm', 'algorithm', srhdSimulation.methods);	//this will change as 'simulations' changes
buildSelect('ADM_algorithm', 'algorithm', admSimulation.methods);	//this will change as 'simulations' changes

buildSelect('eigenDecomposition', 'eigenDecomposition', {Analytic:true, Numeric:true});

buildSelect('explicitMethod', 'explicitMethod', explicitMethods);

{
	const id = 'simulation';
	const key = 'simulation';
	const map = simulations;
	const select = ids[id];
	for (let k in map) {
		const option = Option({innerText : k, appendTo : select});
		if (hydro.state[key] == k) {
			option.setAttribute('selected', 'true');
		}
	}
	select.addEventListener('change', e => {
		const val = select.value;
		hydro.state[key] = val; 
		
		['ADM_controls', 'SRHD_controls', 'Euler_controls'].forEach(id => { hide(ids[id]); });
		show(ids[val+'_controls']);
	
		//...and now select the visible select's selection
		hydro.state.algorithm = ids[val+'_algorithm'].value;
	});
	
	select.value = hydro.state.simulation;
	select.dispatchEvent(new Event('change'));
}

ids.timeStepCFLBased.addEventListener('change', e => {
	if (!ids.timeStepCFLBased.checked) return;
	useCFL = true;
});
ids.timeStepCFL.value = hydro.state.cfl;
ids.timeStepCFL.addEventListener('change', e => {
	const v = Number(ids.timeStepCFL.value);
	if (!isFinite(v)) return;
	hydro.state.cfl = v;
});
ids.timeStepFixed.addEventListener('change', e => {
	if (!ids.timeStepFixed.checked) return;
	useCFL = false;
});
ids.timeStepValue.value = fixedDT;
ids.timeStepValue.addEventListener('change', e => {
	const v = Number(ids.timeStepValue);
	if (!isFinite(v)) return;	//stupid javascript ... convert anything it doesn't understand to NaNs...
	fixedDT = v;
});

const axisColor = [.75, .75, .75, 1];
(new glutil.SceneObject({
	mode : gl.LINES,
	attrs : {
		vertex : new glutil.ArrayBuffer({
			dim : 2, 
			data : [xmin, 0, xmax, 0]
		})
	},
	shader : plainShader,
	uniforms : {
		color : axisColor
	}
})).prependTo(glutil.scene.root);

(new glutil.SceneObject({
	mode : gl.LINES,
	attrs : {
		vertex : new glutil.ArrayBuffer({
			dim : 2, 
			data : [0, ymin, 0, ymax]
		})
	},
	shader : plainShader,
	uniforms : {
		color : axisColor
	}
})).prependTo(glutil.scene.root);


//make static grid
//...might not be so static...
const grid = [];
for (let i = Math.floor(xmin / gridstep) * gridstep; 
	i <= Math.ceil(xmax / gridstep) * gridstep; 
	i += gridstep) 
{
	grid.push(i);
	grid.push(ymin);
	grid.push(i);
	grid.push(ymax);
}
for (let j = Math.floor(ymin / gridstep) * gridstep; 
	j < Math.ceil(ymax / gridstep) * gridstep; 
	j += gridstep)
{
	grid.push(xmin);
	grid.push(j);
	grid.push(xmax);
	grid.push(j);
}
gridObj = new glutil.SceneObject({
	mode : gl.LINES,
	attrs : {
		vertex : new glutil.ArrayBuffer({
			data : grid, 
			dim : 2
		})
	},
	shader : plainShader,
	uniforms : {
		color : [.25,.25,.25,1]
	}
});
gridObj.prependTo(glutil.scene.root);



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
onresize();
window.addEventListener('resize', onresize);
update();
