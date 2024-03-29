function makeExplicitMethods(args) {
const copyState = args.copyState;
const addMulState = args.addMulState;
//depends on copyState(src, dest), addMulState(state, diffState, dt)
const explicitMethods = {
	Euler : function(dt, deriv) {
		let dq_dt = this.tmpq[0];
		copyState(this.q, dq_dt);
		deriv.call(this, dt, dq_dt);
		addMulState(this.q, dq_dt, dt);
	},
	RK2 : function(dt, deriv) {
		let src = this.tmpq[0];
		copyState(this.q, src);
		let k1 = this.tmpq[1];
		copyState(this.q, k1);
		deriv.call(this, dt, k1);
		addMulState(this.q, k1, .5 * dt);
		let k2 = this.tmpq[2];
		copyState(this.q, k2);
		deriv.call(this, dt, k2);
		copyState(src, this.q);
		addMulState(this.q, k2, dt);
	},
	RK4 : function(dt, deriv) {
		let src = this.tmpq[0];
		copyState(this.q, src);
		let k1 = this.tmpq[1];
		copyState(this.q, k1);
		deriv.call(this, dt, k1);
		addMulState(this.q, k1, .5 * dt);
		let k2 = this.tmpq[2];
		copyState(this.q, k2);
		deriv.call(this, dt, k2);
		copyState(src, this.q);
		addMulState(this.q, k2, .5 * dt);
		let k3 = this.tmpq[3];
		copyState(this.q, k3);
		deriv.call(this, dt, k3);
		copyState(src, this.q);
		addMulState(this.q, k3, dt);
		let k4 = this.tmpq[4];
		copyState(this.q, k4);
		deriv.call(this, dt, k4);
		copyState(src, this.q);
		addMulState(this.q, k1, dt / 6);
		addMulState(this.q, k2, dt / 3);
		addMulState(this.q, k3, dt / 3);
		addMulState(this.q, k4, dt / 6);
	},
	ICN3 : function(dt, deriv) {
		//first iteration
		let srcQ = this.tmpq[0];
		copyState(this.q, srcQ);
		let firstK = this.tmpq[1];
		copyState(this.q, firstK);
		deriv.call(this, dt, firstK);
		addMulState(this.q, firstK, dt);
		let k = this.tmpq[2];
		copyState(this.q, k);

		//second and so on
		for (let i = 1; i < 3; ++i) {
			deriv.call(this, dt, k);
			copyState(srcQ, this.q);
			addMulState(this.q, k, .5 * dt);
			addMulState(this.q, firstK, .5 * dt);
		}
	}
};
return explicitMethods;
}
export {makeExplicitMethods};
