#![allow(dead_code, unused_macros)]

use std::{fs::File, io::Write};

#[derive(Clone)]
struct Field {
    u: Box<[f64]>,
    v: Box<[f64]>,
    p: Box<[f64]>,
    e: Box<[f64]>,
    rho: Box<[f64]>,
    k: Box<[f64]>,
    cp: Box<[f64]>,
}

impl Field {
    fn new_zeros(size: usize) -> Self {
        let b = vec![0.0; size].into_boxed_slice();
        Self{
            u: b.clone(),
            v: b.clone(),
            p: b.clone(),
            e: b.clone(),
            rho: b.clone(),
            k: b.clone(),
            cp: b,
        }
    }

    fn set_initcond(&mut self, idx: usize, ip: &InitParams) {
        self.u[idx] = 0.0;
        self.v[idx] = 0.0; 
        self.p[idx] = ip.p; 
        self.e[idx] = ip.e; 
        self.rho[idx] = ip.rho; 
        self.k[idx] = ip.k;
        self.cp[idx] = ip.cp;
    }

    fn temperature_at(&self, idx: usize) -> f64 {
        let k = self.k[idx];
        (self.p[idx] * k) / (self.rho[idx] * self.cp[idx] * (k - 1.0))
    }
}

#[derive(Clone)]
struct BndField {
    x: Box<[f64]>,
    p: Box<[f64]>,

    rho_x: Box<[f64]>,
    rho_y_x: Box<[f64]>,
    rho_x_x: Box<[f64]>,
    rho_e_x: Box<[f64]>,
    rho_k_x: Box<[f64]>,
    rho_cp_x: Box<[f64]>,
}

impl BndField {
    fn new_zeros(size: usize) -> Self {
        let b = vec![0.0; size].into_boxed_slice();
        Self{
            x: b.clone(),
            p: b.clone(),
            rho_x: b.clone(),
            rho_y_x: b.clone(),
            rho_x_x: b.clone(),
            rho_e_x: b.clone(),
            rho_k_x: b.clone(),
            rho_cp_x: b,
        }
    }
}

#[derive(Clone)]
struct Solver {
    field: Field,
    bnds_z: BndField,
    bnds_r: BndField,

    dtime: f64,
    dz: f64,
    dr: f64,

    num_z: usize,
    num_r: usize, 
}

macro_rules! update_avg_bnds_along_z {
    ($self:ident, $bnd:ident, $v0:ident, $f:ident, $v1:ident) => {{
        for j in 0..$self.num_r {
            let mut bidx = j * ($self.num_z + 1);
            let mut fidx = (j + 1) * ($self.num_z + 2);
            for _ in 0..$self.num_z + 1 {
                $self.$bnd.$v0[bidx] = 0.5 * ($self.$f.$v1[fidx] + $self.$f.$v1[fidx + 1]);
                bidx += 1;
                fidx += 1;
            }
        }
    }};
    ($self:ident, $idx:ident, $idxr:ident, $br:ident, $expression:expr) => {{
        for j in 1..$self.num_r + 1 {
            let mut $br = (j - 1) * ($self.num_z + 1);
            let mut $idx = j * ($self.num_z + 2);
            for _ in 0..$self.num_z + 1 {
                let $idxr = $idx + 1;
        
                $expression
        
                $br += 1;
                $idx += 1;
            }
        }
    }};
}

macro_rules! update_avg_bnds_along_r {
    ($self:ident, $bnd:ident, $v0:ident, $f:ident, $v1:ident) => {{
        for j in 0..$self.num_r + 1 {
            let mut bidx = j * $self.num_z;
            let mut fidx = j * ($self.num_z + 2) + 1;
            let dfv = $self.num_z + 2;
            for _ in 0..$self.num_z {
                $self.$bnd.$v0[bidx] = 0.5 * ($self.$f.$v1[fidx] + $self.$f.$v1[fidx + dfv]);
                bidx += 1;
                fidx += 1;
            }
        }
    }};
    ($self:ident, $idx:ident, $idxt:ident, $bt:ident, $expression:expr) => {{
        for j in 0..$self.num_r + 1 {
            let mut $bt = j * $self.num_z;
            let mut $idx = j * ($self.num_z + 2) + 1;
            for _ in 0..$self.num_z {
                let $idxt = $idx + $self.num_z + 2;
        
                $expression
        
                $bt += 1;
                $idx += 1;
            }
        }
    }};
}

macro_rules! approx_loop {
    ($self:ident, $j:ident, $idx:ident,
        $bl:ident, $br:ident, $bb:ident, $bt:ident, 
        $expression:expr) => {{

        for $j in 1..$self.num_r + 1 {
            let mut $bl = ($j - 1) * ($self.num_z + 1);
            let mut $bb = ($j - 1) * $self.num_z;
            let mut $idx = ($self.num_z + 2) * $j + 1;
            for _ in 1..$self.num_z + 1 {
                let $br = $bl + 1;
                let $bt = $bb + $self.num_z;
        
                $expression
        
                $bl += 1;
                $bb += 1;
                $idx += 1;
            }
        }
    }};
}

macro_rules! update_bnd_prod_along_z {
    ($self:ident, $idx:ident, $br:ident) => {{
        let rho_bur = $self.field.rho[$idx] * $self.bnds_z.x[$br];
        let upd = |prod: &mut [f64], val: &[f64]| prod[$br] = rho_bur * val[$idx];
        $self.bnds_z.rho_x[$br] = rho_bur;
        upd(&mut $self.bnds_z.rho_x_x, &$self.field.u);
        upd(&mut $self.bnds_z.rho_y_x, &$self.field.v);
        upd(&mut $self.bnds_z.rho_e_x, &$self.field.e);
        upd(&mut $self.bnds_z.rho_k_x, &$self.field.k);
        upd(&mut $self.bnds_z.rho_cp_x, &$self.field.cp);
    }};
}

macro_rules! update_bnd_prod_along_r {
    ($self:ident, $idx:ident, $br:ident) => {{
        let rho_bur = $self.field.rho[$idx] * $self.bnds_r.x[$br];
        let upd = |prod: &mut [f64], val: &[f64]| prod[$br] = rho_bur * val[$idx];
        $self.bnds_r.rho_x[$br] = rho_bur;
        upd(&mut $self.bnds_r.rho_x_x, &$self.field.v);
        upd(&mut $self.bnds_r.rho_y_x, &$self.field.u);
        upd(&mut $self.bnds_r.rho_e_x, &$self.field.e);
        upd(&mut $self.bnds_r.rho_k_x, &$self.field.k);
        upd(&mut $self.bnds_r.rho_cp_x, &$self.field.cp);
    }};
}

macro_rules! div_rho_x_w {
    ($self:ident, $j:ident, 
        $bl:ident, $br:ident, $bb:ident, $bt:ident, 
        $dt_div_dz:ident, $along_z:ident, 
        $dt_div_sdr:ident, $along_r:ident) => {{
        
        $dt_div_dz * ($self.bnds_z.$along_z[$br] - $self.bnds_z.$along_z[$bl])
        - $dt_div_sdr * ($j as f64 * $self.bnds_r.$along_r[$bt] - ($j - 1) as f64 * $self.bnds_r.$along_r[$bb])
    }};
}

macro_rules! update_bnd_conds_scalar {
    ($self:ident, $scalar:ident) => {{
        for i in 1..$self.num_z + 1 {
            $self.field.$scalar[i] = $self.field.$scalar[$self.num_z + 2 + i];
        }
        let beg_top = ($self.num_r + 1) * ($self.num_z + 2);
        let beg_under_top = beg_top - ($self.num_z + 2);
        for i in 1..$self.num_z + 1 {
            $self.field.$scalar[beg_top + i] = $self.field.$scalar[beg_under_top + i];
        }

        for j in 1..$self.num_r + 1 {
            let beg = j * ($self.num_z + 2);
            $self.field.$scalar[beg] = $self.field.$scalar[beg + 1];
            $self.field.$scalar[beg + $self.num_z + 1] = $self.field.$scalar[beg + $self.num_z];
        }        
    }};
}

macro_rules! update_bnd_conds_velocity {
    ($self:ident) => {{
        for i in 1..$self.num_z + 1 {
            $self.field.u[i] = $self.field.u[$self.num_z + 2 + i];
            $self.field.v[i] = -$self.field.v[$self.num_z + 2 + i];
        }
        let beg_top = ($self.num_r + 1) * ($self.num_z + 2);
        let beg_under_top = beg_top - ($self.num_z + 2);
        for i in 1..$self.num_z + 1 {
            $self.field.u[beg_top + i] = $self.field.u[beg_under_top + i];
            $self.field.v[beg_top + i] = -$self.field.v[beg_under_top + i];
        }

        for j in 1..$self.num_r + 1 {
            let beg = j * ($self.num_z + 2);
            $self.field.u[beg] = -$self.field.u[beg + 1];
            $self.field.u[beg + $self.num_z + 1] = -$self.field.u[beg + $self.num_z];
            $self.field.v[beg] = $self.field.v[beg + 1];
            $self.field.v[beg + $self.num_z + 1] = $self.field.v[beg + $self.num_z];
        }        
    }};
}

impl Solver {
    fn init_conds(
        air: &InitParams, fuel: &InitParams, 
        num_r: usize, num_z_hp: usize, num_z: usize
    ) -> Field {
        let mut field = Field::new_zeros((num_z + 2) * (num_r + 2));
        for j in 0..num_r + 2 {
            let hbeg = j * (num_z + 2);
            let hend = hbeg + num_z_hp + 1;
            let lend = hbeg + num_z + 2;

            field.u[hbeg..lend].fill(0.0);
            field.v[hbeg..lend].fill(0.0);
            field.p[hbeg..hend].fill(fuel.p);
            field.p[hend..lend].fill(air.p);
            field.e[hbeg..hend].fill(fuel.e);
            field.e[hend..lend].fill(air.e);
            field.rho[hbeg..hend].fill(fuel.rho);
            field.rho[hend..lend].fill(air.rho);
            field.k[hbeg..hend].fill(fuel.k);
            field.k[hend..lend].fill(air.k);
            field.cp[hbeg..hend].fill(fuel.cp);
            field.cp[hend..lend].fill(air.cp);
        }
        field
    }

    fn new(
        air: &InitParams, fuel: &InitParams, 
        num_r: usize, num_z_hp: usize, num_z_lp: usize, 
        dz: f64, dr: f64
    ) -> Self {

        let num_z = num_z_hp + num_z_lp;
        let field = Self::init_conds(air, fuel, num_r, num_z_hp, num_z);
        let bnds_z = BndField::new_zeros((num_z + 1) * num_r);
        let bnds_r = BndField::new_zeros(num_z * (num_r + 1));

        Solver{ 
            field, bnds_z, bnds_r, 
            dtime: 0.0, dz, dr, 
            num_r, num_z
        }
    }

    fn update_bnds_p(&mut self) {
        update_avg_bnds_along_z!(self, bnds_z, p, field, p);
        update_avg_bnds_along_r!(self, bnds_r, p, field, p);
    }

    fn update_bnds_uv(&mut self) {
        update_avg_bnds_along_z!(self, bnds_z, x, field, u);
        update_avg_bnds_along_r!(self, bnds_r, x, field, v);
    }

    fn update_pressure(&mut self) {
        for j in 1..self.num_r + 1 {
            let mut idx = (self.num_z + 2) * j + 1;
            for _ in 1..self.num_z + 1 {
                self.field.p[idx] = 
                    (self.field.k[idx] - 1.0) * self.field.rho[idx] 
                    * (self.field.e[idx] - 0.5 
                    * (self.field.u[idx].powi(2) + self.field.v[idx].powi(2)));
                idx += 1;
            }
        }
    }

    fn euler(&mut self, dtime: f64) -> &mut Self {
        self.dtime = dtime;
        self.update_bnds_p();
        approx_loop!(self, j, idx,
            bl, br, bb, bt, {
                let pl = self.bnds_z.p[bl];
                let pr = self.bnds_z.p[br];
                let pb = self.bnds_r.p[bb];
                let pt = self.bnds_r.p[bt];

                let dt_div_rho = self.dtime / self.field.rho[idx];
                self.field.u[idx] -= (pr - pl) / self.dz * dt_div_rho;
                self.field.v[idx] -= (pt - pb) / self.dr * dt_div_rho;
            }
        );
        update_bnd_conds_velocity!(self);
        self.update_bnds_uv();

        approx_loop!(self, j, idx,
            bl, br, bb, bt, {
                let pl = self.bnds_z.p[bl];
                let pr = self.bnds_z.p[br];
                let pb = self.bnds_r.p[bb];
                let pt = self.bnds_r.p[bt];

                let ul = self.bnds_z.x[bl];
                let ur = self.bnds_z.x[br];
                let vb = self.bnds_r.x[bb];
                let vt = self.bnds_r.x[bt];

                let dt_div_rho = self.dtime / self.field.rho[idx];
                self.field.e[idx] -= (
                    (pr * ur - pl * ul) / self.dz 
                    + (j as f64 * pt * vt - (j - 1) as f64 * pb * vb
                    ) / ((j as f64 - 0.5) * self.dr)
                ) * dt_div_rho;
            }
        );
        update_bnd_conds_scalar!(self, e);

        self
    }

    fn lagrange(&mut self) -> &mut Self {
        update_avg_bnds_along_z!(self, idx, idxr, br, {
                if self.bnds_z.x[br] > 0.0 {
                    update_bnd_prod_along_z!(self, idx, br);
                } else {
                    update_bnd_prod_along_z!(self, idxr, br);
                }
            }
        );
        update_avg_bnds_along_r!(self, idx, idxt, bt, {
                if self.bnds_r.x[bt] > 0.0 {
                    update_bnd_prod_along_r!(self, idx, bt);
                } else {
                    update_bnd_prod_along_r!(self, idxt, bt);
                }
            }
        );
        self
    }

    fn finalize(&mut self) -> &mut Self {
        approx_loop!(self, j, idx,
            bl, br, bb, bt, {
                let rho_euler = self.field.rho[idx];
                let dt_div_dz = self.dtime / self.dz;
                let dt_div_sdr = self.dtime / ((j as f64 - 0.5) * self.dr);

                self.field.rho[idx] -= 
                    div_rho_x_w!(self, j, bl, br, bb, bt, 
                        dt_div_dz, rho_x, dt_div_sdr, rho_x);

                let inv_rho = 1.0 / self.field.rho[idx];

                self.field.u[idx] = (
                    rho_euler * self.field.u[idx] 
                    - div_rho_x_w!(self, j, bl, br, bb, bt, 
                        dt_div_dz, rho_x_x, dt_div_sdr, rho_y_x)
                ) * inv_rho;

                self.field.v[idx] = (
                    rho_euler * self.field.v[idx] 
                    - div_rho_x_w!(self, j, bl, br, bb, bt, 
                        dt_div_dz, rho_y_x, dt_div_sdr, rho_x_x)
                ) * inv_rho;

                self.field.e[idx] = (
                    rho_euler * self.field.e[idx] 
                    - div_rho_x_w!(self, j, bl, br, bb, bt, 
                        dt_div_dz, rho_e_x, dt_div_sdr, rho_e_x)
                ) * inv_rho;

                self.field.k[idx] = (
                    rho_euler * self.field.k[idx] 
                    - div_rho_x_w!(self, j, bl, br, bb, bt, 
                        dt_div_dz, rho_k_x, dt_div_sdr, rho_k_x)
                ) * inv_rho;

                self.field.cp[idx] = (
                    rho_euler * self.field.cp[idx] 
                    - div_rho_x_w!(self, j, bl, br, bb, bt, 
                        dt_div_dz, rho_cp_x, dt_div_sdr, rho_cp_x)
                ) * inv_rho;
            }
        );
        self.update_pressure();
        update_bnd_conds_velocity!(self);
        update_bnd_conds_scalar!(self, rho);
        update_bnd_conds_scalar!(self, e);
        update_bnd_conds_scalar!(self, k);
        update_bnd_conds_scalar!(self, cp);
        update_bnd_conds_scalar!(self, p);

        self
    }

    fn output_by<T>(
        &self, out: &mut impl Write, fvals: &[T],
        f: impl Fn(&mut dyn Iterator<Item=&T>) -> String
    ) {
        for j in 1..1 + self.num_r {
            let beg = (self.num_z + 2) * j + 1;
            writeln!(
                out, "{}", 
                f(&mut fvals.iter().skip(beg).take(self.num_z))
            ).unwrap();
        }
    }

    fn output_pressure(&self, out: &mut impl Write) {
        self.output_by(out, &self.field.p, |iter| {
            iter.map(|p| format!("{:e}", p))
                .reduce(|acc, p| format!("{} {}", acc, p))
                .unwrap()
        })
    }

    fn output_density(&self, out: &mut impl Write) {
        self.output_by(out, &self.field.rho, |iter| {
            iter.map(|rho| format!("{:e}", rho))
                .reduce(|acc, rho| format!("{} {}", acc, rho))
                .unwrap()
        })
    }

    fn output_heat_capacity_ratio(&self, out: &mut impl Write) {
        self.output_by(out, &self.field.k, |iter| {
            iter.map(|k| format!("{:e}", k))
                .reduce(|acc, k| format!("{} {}", acc, k))
                .unwrap()
        })
    }

    fn output_temperature(&self, out: &mut impl Write) {
        let ts: Vec<f64> = (0..self.field.p.len()).map(|x| self.field.temperature_at(x)).collect();
        self.output_by(out, &ts, |iter| {
            iter.map(|t| format!("{:e}", t))
                .reduce(|acc, t| format!("{} {}", acc, t))
                .unwrap()
        })
    }

    fn output_velocity(&self, out: &mut impl Write) {
        let uv: Vec<_> = self.field.u.iter().zip(self.field.v.iter()).collect();
        self.output_by(out, &uv, |iter| {
            iter.map(|x| format!("{:e} {:e}", x.0, x.1))
                .reduce(|acc, x| format!("{} {}", acc, x))
                .unwrap()
        })
    }

    fn output_everything(&self, pout: &mut impl Write, dout: &mut impl Write, kout: &mut impl Write, tout: &mut impl Write, vout: &mut impl Write) {
        self.output_pressure(pout);
        self.output_density(dout);
        self.output_heat_capacity_ratio(kout);
        self.output_temperature(tout);
        self.output_velocity(vout);
    }
}

const ATM: f64 = 101325.0;

#[derive(Clone, Copy, Debug)]
struct InitParams {
    p: f64,
    t: f64,
    r: f64,
    k: f64,
    cp: f64,
    rho: f64,
    e: f64,
}

impl InitParams {
    fn new(p: f64, t: f64, r: f64, k: f64) -> Self {
        let cp = r * k / (k - 1.0);
        let rho = p / (r * t);
        let e = p / ((k - 1.0) * rho);
        InitParams{ p, t, r, k, cp, rho, e }
    }
}

fn air_init_params() -> InitParams {
    InitParams::new(ATM, 293.0, 287.0, 1.4)
}

fn fuel_init_params() -> InitParams {
    InitParams::new(40.0 * ATM, 400.0, 310.0, 1.2)
}

fn main() {
    let air = air_init_params();
    let fuel = fuel_init_params();
    println!("Air {:?}", air);
    println!("Fuel {:?}", fuel);
    let num_r = 10;
    let num_z_hp = 50;
    let num_z_lp = 450;
    let dtime = 1e-6;
    let dz = 5e-3;
    let dr = dz;

    let mut solver = Solver::new(&air, &fuel, num_r, num_z_hp, num_z_lp, dz, dr);

    let mut pf = File::create("pressure.ssv").unwrap();
    let mut df = File::create("density.ssv").unwrap();
    let mut kf = File::create("hc-ratio.ssv").unwrap();
    let mut tf = File::create("temperature.ssv").unwrap();
    let mut vf = File::create("velocity.ssv").unwrap();
    solver.output_everything(&mut pf, &mut df, &mut kf, &mut tf, &mut vf);
    for i in 0..20000 {
        solver.euler(dtime).lagrange().finalize();
        if (i + 1) % 200 == 0 {
            [&mut pf, &mut df, &mut kf, &mut tf, &mut vf].iter_mut()
                .for_each(|f| f.write_fmt(format_args!("\n")).unwrap());
            solver.output_everything(&mut pf, &mut df, &mut kf, &mut tf, &mut vf);
        }
    }
}
