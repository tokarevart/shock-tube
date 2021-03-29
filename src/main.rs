#![allow(dead_code)]

use std::{fs::File, io::Write};
use itertools::izip;

// #[derive(Clone, Copy, Debug)]
// struct FieldPoint {
//     u: f64,
//     v: f64,
//     p: f64,
//     e: f64,
//     rho: f64,
//     k: f64,
//     cp: f64,
// }

// impl FieldPoint {
//     fn zeros() -> Self {
//         FieldPoint{ 
//             u: 0.0, v: 0.0, 
//             p: 0.0, e: 0.0, rho: 0.0, 
//             k: 0.0, cp: 0.0 
//         }
//     }

//     fn from_initconds(ip: &InitParams) -> Self {
//         FieldPoint{ 
//             u: 0.0, v: 0.0, 
//             p: ip.p, e: ip.e, rho: ip.rho, 
//             k: ip.k, cp: ip.cp 
//         }
//     }

//     fn temperature(&self) -> f64 {
//         (self.p * self.k) / (self.rho * self.cp * (self.k - 1.0))
//     }
// }

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
    // fn iter(&self) -> impl Iterator<Item=
    //         (&f64, &f64, &f64, &f64, &f64, &f64, &f64)> {
    //     izip!(
    //         self.u.iter(), 
    //         self.v.iter(), 
    //         self.p.iter(), 
    //         self.e.iter(), 
    //         self.rho.iter(), 
    //         self.k.iter(), 
    //         self.cp.iter(), 
    //     )
    // }

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

    time: f64,
    dtime: f64,
    dz: f64,
    dr: f64,

    num_z: usize,
    num_r: usize, 
}

macro_rules! update_bnds_along_z {
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
}

macro_rules! update_bnds_along_r {
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
        dtime: f64, dz: f64, dr: f64
    ) -> Self {

        let num_z = num_z_hp + num_z_lp;
        let field = Self::init_conds(air, fuel, num_r, num_z_hp, num_z);
        let bnds_z = BndField::new_zeros((num_z + 1) * num_r);
        let bnds_r = BndField::new_zeros(num_z * (num_r + 1));

        Solver{ 
            field, bnds_z, bnds_r, 
            time: 0.0, dtime, dz, dr, 
            num_r, num_z
        }
    }

    fn update_bnds_p(&mut self) {
        // for j in 0..self.num_r {
        //     let mut bidx = j * (self.num_z + 1);
        //     let mut fidx = (j + 1) * (self.num_z + 2);
        //     for i in 0..self.num_z + 1 {
        //         self.bnds_z.p[bidx] = 0.5 * (self.field.p[fidx] + self.field.p[fidx + 1]);
        //         bidx += 1;
        //         fidx += 1;
        //     }
        // }
        update_bnds_along_z!(self, bnds_z, p, field, p);

        // for j in 0..self.num_r + 1 {
        //     let mut bidx = j * self.num_z;
        //     let mut fidx = j * (self.num_z + 2) + 1;
        //     let dfv = self.num_z + 2;
        //     for i in 0..self.num_z {
        //         self.bnds_r.p[bidx] = 0.5 * (self.field.p[fidx] + self.field.p[fidx + dfv]);
        //         bidx += 1;
        //         fidx += 1;
        //     }
        // }
        update_bnds_along_r!(self, bnds_r, p, field, p);
    }

    fn update_bnds_uv(&mut self) {
        // for j in 0..self.num_r {
        //     let mut bidx = j * (self.num_z + 1);
        //     let mut fidx = (j + 1) * (self.num_z + 2);
        //     for i in 0..self.num_z + 1 {
        //         self.bnds_z.x[bidx] = 0.5 * (self.field.u[fidx] + self.field.u[fidx + 1]);
        //         bidx += 1;
        //         fidx += 1;
        //     }
        // }
        update_bnds_along_z!(self, bnds_z, x, field, u);

        // for j in 0..self.num_r + 1 {
        //     let mut bidx = j * self.num_z;
        //     let mut fidx = j * (self.num_z + 2) + 1;
        //     let dfv = self.num_z + 2;
        //     for i in 0..self.num_z {
        //         self.bnds_r.x[bidx] = 0.5 * (self.field.v[fidx] + self.field.v[fidx + dfv]);
        //         bidx += 1;
        //         fidx += 1;
        //     }
        // }
        update_bnds_along_r!(self, bnds_r, x, field, v);
    }

    fn euler(&mut self) {
        let num_r = self.num_r;
        let num_z = self.num_z;

        self.update_bnds_p();

        for j in 1..num_r + 1 {
            let mut bidxl = (j - 1) * (num_z + 1) - 1;
            let mut bidxb = (j - 1) * num_z - 1;
            let mut idx = (num_z + 2) * j + 1;
            for i in 1..num_z + 1 {
                let bidxr = bidxl + 1;
                let bidxt = bidxb + num_z;
                let pl = self.bnds_z.p[bidxl];
                let pr = self.bnds_z.p[bidxr];
                let pb = self.bnds_r.p[bidxb];
                let pt = self.bnds_r.p[bidxt];
                let dt_div_rho = self.dtime / self.field.rho[idx];
                self.field.u[idx] -= (pr - pl) / self.dz * dt_div_rho;
                self.field.v[idx] -= (pt - pb) / self.dr * dt_div_rho;
                self.update_bnds_uv();
                // self.field.e[idx] -= () * dt_div_rho;

                bidxl += 1;
                bidxb += 1;
                idx += 1;
            }
        }
        //
    }

    fn output_by<T>(
        &self, out: &mut impl Write, fvals: &[T],
        f: impl Fn(&mut dyn Iterator<Item=&T>) -> String
    ) {

        for j in 1..1 + self.num_r {
            let beg = (self.num_z + 2) * j + 1;
            let end = beg + self.num_z;
            writeln!(
                out, "{}", 
                f(&mut fvals.iter().skip(beg).take(self.num_z))
            ).unwrap();
        }
    }

    fn output_pressure(&self, out: &mut impl Write) {
        self.output_by(out, &self.field.p, |iter| {
            iter.map(|p| format!("{:e}", p))
                .reduce(|acc, p| format!("{},{}", acc, p))
                .unwrap()
        })
    }

    fn output_temperature(&self, out: &mut impl Write) {
        let ts: Vec<f64> = (0..self.field.p.len()).map(|x| self.field.temperature_at(x)).collect();
        self.output_by(out, &ts, |iter| {
            iter.map(|t| format!("{:e}", t))
                .reduce(|acc, t| format!("{},{}", acc, t))
                .unwrap()
        })
    }

    fn output_velocity(&self, out: &mut impl Write) {
        let uv: Vec<_> = self.field.u.iter().zip(self.field.v.iter()).collect();
        self.output_by(out, &uv, |iter| {
            iter.map(|x| format!("{:e}|{:e}", x.0, x.1))
                .reduce(|acc, x| format!("{},{}", acc, x))
                .unwrap()
        })
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
    let num_r = 10;
    let num_z_hp = 50;
    let num_z_lp = 450;
    let dtime = 1e-6;
    let dz = 5e-3;
    let dr = dz;

    let solver = Solver::new(&air, &fuel, num_r, num_z_hp, num_z_lp, dtime, dz, dr);

    solver.output_pressure(&mut File::create("pressure.csv").unwrap());
    solver.output_temperature(&mut File::create("temperature.csv").unwrap());
    solver.output_velocity(&mut File::create("velocity.csv").unwrap());

    println!("{:?}", air);
    println!("{:?}", fuel);
}
