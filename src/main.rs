#![allow(dead_code)]

use std::{fs::File, io::Write};

#[derive(Clone, Copy, Debug)]
struct FieldPoint {
    u: f64,
    v: f64,
    p: f64,
    e: f64,
    rho: f64,
    k: f64,
    cp: f64,
}

impl FieldPoint {
    fn zeros() -> Self {
        FieldPoint{ 
            u: 0.0, v: 0.0, 
            p: 0.0, e: 0.0, rho: 0.0, 
            k: 0.0, cp: 0.0 
        }
    }

    fn from_initconds(ip: &InitParams) -> Self {
        FieldPoint{ 
            u: 0.0, v: 0.0, 
            p: ip.p, e: ip.e, rho: ip.rho, 
            k: ip.k, cp: ip.cp 
        }
    }

    fn temperature(&self) -> f64 {
        (self.p * self.k) / (self.rho * self.cp * (self.k - 1.0))
    }
}

struct Solver {
    field: Box<[FieldPoint]>,
    lagr_bnds_u: Box<[FieldPoint]>,
    lagr_bnds_v: Box<[FieldPoint]>,

    time: f64,
    dtime: f64,
    dz: f64,
    dr: f64,

    num_z: usize,
    num_r: usize, 
}

impl Solver {
    fn field_get(&self, i: usize, j: usize) -> &FieldPoint {
        &self.field[j * self.num_z + i]
    }

    fn init_conds(
        air: &InitParams, fuel: &InitParams, 
        num_r: usize, num_z_hp: usize, num_z: usize
    ) -> Box<[FieldPoint]> {

        let mut field = vec![
            FieldPoint::from_initconds(air); 
            (num_z + 2) * (num_r + 2)
        ].into_boxed_slice();
        let fuel_fp = FieldPoint::from_initconds(fuel);
        for j in 0..num_r + 2 {
            let beg = j * (num_z + 2);
            let end = beg + num_z_hp + 1;
            field[beg..end].fill(fuel_fp);
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
        let lagr_bnds_u = vec![
            FieldPoint::zeros(); 
            (num_z + 1) * num_r
        ].into_boxed_slice();
        let lagr_bnds_v = lagr_bnds_u.clone();

        Solver{ 
            field, lagr_bnds_u, lagr_bnds_v, 
            time: 0.0, dtime, dz, dr, 
            num_r, num_z
        }
    }

    fn output_by(
        &self, out: &mut impl Write, 
        f: impl Fn(std::slice::Iter<'_, FieldPoint>
    ) -> String) {

        for j in 1..1 + self.num_r {
            let beg = (self.num_z + 2) * j + 1;
            let end = beg + self.num_z;
            writeln!(
                out, "{}", 
                f(self.field[beg..end].iter())
            ).unwrap();
        }
    }

    fn output_pressure(&self, out: &mut impl Write) {
        self.output_by(out, |iter| {
            iter.map(|x| format!("{:e}", x.p))
                .reduce(|acc, p| format!("{},{}", acc, p))
                .unwrap()
        })
    }

    fn output_temperature(&self, out: &mut impl Write) {
        self.output_by(out, |iter| {
            iter.map(|x| format!("{:e}", x.temperature()))
                .reduce(|acc, p| format!("{},{}", acc, p))
                .unwrap()
        })
    }

    fn output_velocity(&self, out: &mut impl Write) {
        self.output_by(out, |iter| {
            iter.map(|x| format!("{:e}^{:e}", x.u, x.v))
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
