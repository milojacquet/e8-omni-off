use crate::e8::E8;
use crate::point::D8;
use crate::point::Orbit;
use crate::point::Point;
use fxhash::FxHashMap;
use fxhash::FxHashSet;

use crate::Mirror;
use crate::e8::MirrorSet;
use std::collections::VecDeque;
use std::fmt::Display;
use std::io::Write;

fn num_length_i16(x: i16) -> u64 {
    x.to_string().len() as u64
}

fn num_length_usize(x: usize) -> u64 {
    x.to_string().len() as u64
}

fn num_length_u64(x: u64) -> u64 {
    x.to_string().len() as u64
}

#[derive(Debug, Clone)]
struct PointSet {
    orbits: Vec<(Orbit, (E8, D8))>,
    lookup: FxHashMap<Orbit, u64>,
    len: u64,
}

// TODO: make the iterator return what E8 you need to get there
impl PointSet {
    fn new(iter: impl Iterator<Item = (Point, E8)>) -> Self {
        let mut orbits = Vec::new();
        let mut lookup = FxHashMap::from_iter([]);
        let mut offset = 0;
        for (point, e8) in iter {
            // point = original * e8
            // point = orbitrep * d8
            // orbitrep = point * d8.inv() = original * e8 * d8.inv()
            orbits.push((point.orbit, (e8, point.d8.inv())));
            lookup.insert(point.orbit, offset);
            offset += point.orbit.size();
        }
        Self {
            orbits,
            lookup,
            len: offset,
        }
    }

    fn len(&self) -> u64 {
        self.len
    }

    fn index(&self, point: Point) -> u64 {
        let offset = self.lookup[&point.orbit];
        offset + point.orbit_index()
    }

    fn iter(&self) -> impl Iterator<Item = (Point, (E8, D8))> {
        self.orbits
            .iter()
            .flat_map(|(orbit, (e8, d8))| orbit.iter().map(|point| (point, (*e8, *d8 * point.d8))))
    }
}

pub fn write_spaced<T: Display>(
    mut writer: impl Write,
    iter: impl Iterator<Item = T>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut is_first = true;
    for item in iter {
        if is_first {
            is_first = false
        } else {
            write!(writer, " ")?;
        }
        write!(writer, "{}", item)?;
    }
    Ok(())
}

impl MirrorSet {
    pub fn face_types(self) -> [Vec<Self>; 9] {
        let mut face_types = [const { Vec::new() }; 9];
        for set in Self::iter_all() {
            let components = set.components();
            if components.iter().all(|comp| !(*comp & self).is_empty()) {
                face_types[set.size() as usize].push(set);
            }
        }
        face_types
    }

    pub fn face_center(self, face: Self) -> Self {
        // Add largest possible components to face disjoint from face and self and take complement
        // Equivalent: complement of all nodes in face or (nonadjacent to face and not in self)
        // Equivalent to: all nodes not in face and (adjacent to face or in self)

        Mirror::ALL
            .into_iter()
            .filter(|&m| {
                !face.has_mirror(m)
                    && (self.has_mirror(m)
                        || Mirror::ALL
                            .into_iter()
                            .any(|mm| face.has_mirror(mm) && mm.link(m) == 3))
            })
            .fold(Self::empty(), |x, y| x | Self::from_mirror(y))
    }

    pub fn off_size_estimate(self) -> u64 {
        let mut size = 0;
        let face_types = self.face_types();
        let point_sets: [PointSet; 9] = face_types.clone().map(|types| {
            PointSet::new(
                types
                    .iter()
                    .flat_map(|face| self.face_center(*face).vertex_orbits()),
            )
        });
        // dbg!(face_types, face_centers);

        for (orbit, _) in &point_sets[0].orbits {
            size += (8 * num_length_i16(orbit.rep.max()) + 12) * orbit.size();
        }

        for i in 2..8 {
            for &face_type in face_types[i].iter() {
                // kind of overkill but if it's a problem wait until you see what comes next
                let mut subfaces = FxHashSet::from_iter([]);
                for &subface_type in &face_types[i - 1] {
                    if face_type.contains(subface_type) {
                        let mut stack =
                            VecDeque::from_iter([self.face_center(subface_type).vertex()]);
                        while let Some(v) = stack.pop_front() {
                            if subfaces.insert(v) {
                                for mirror in face_type.mirrors() {
                                    stack.push_back(v * mirror.mat());
                                }
                            }
                        }
                    }
                }

                size += (num_length_usize(subfaces.len())
                    + subfaces.len() as u64 * (1 + num_length_u64(point_sets[i - 1].len())))
                    * self.face_center(face_type).vertex_count();
            }
        }

        size
    }

    pub fn write_off(self, mut writer: impl Write) -> Result<(), Box<dyn std::error::Error>> {
        let face_types = self.face_types();
        let point_sets: [PointSet; 9] = face_types.clone().map(|types| {
            PointSet::new(
                types
                    .iter()
                    .flat_map(|face| self.face_center(*face).vertex_orbits()),
            )
        });
        // dbg!(face_types, face_centers);

        write!(writer, "8OFF\n")?;
        write!(writer, "{} ", point_sets[0].len())?;
        write!(writer, "{} ", point_sets[2].len())?;
        write!(writer, "{} ", point_sets[1].len())?;
        for i in 3..8 {
            // intentionally omitting 8
            write!(writer, "{} ", point_sets[i].len())?;
        }
        write!(writer, "\n\n")?;

        println!("Vertices");
        write!(writer, "# Vertices\n")?;
        for (vertex, _) in point_sets[0].iter() {
            write_spaced(&mut writer, vertex.vec().iter())?;
            write!(writer, "\n")?;
        }
        write!(writer, "\n")?;

        write!(writer, "# Faces\n")?;
        for &face_type in face_types[2].iter() {
            println!("Faces: {face_type:?}");
            let vertex = self.vertex();
            let &[m1, m2] = &face_type.mirrors().collect::<Vec<_>>()[..] else {
                panic!("not two")
            };

            let vertices = if m1.link(m2) == 2 {
                vec![
                    vertex,
                    vertex * m1.mat(),
                    vertex * m2.mat() * m1.mat(),
                    vertex * m2.mat(),
                ]
            } else if self & face_type == face_type {
                vec![
                    vertex,
                    vertex * m1.mat(),
                    vertex * m2.mat() * m1.mat(),
                    vertex * m1.mat() * m2.mat() * m1.mat(),
                    vertex * m1.mat() * m2.mat(),
                    vertex * m2.mat(),
                ]
            } else {
                vec![
                    vertex,
                    vertex * m2.mat() * m1.mat(),
                    vertex * m1.mat() * m2.mat(),
                ]
            };

            // already did this once, is that ok
            for (_point, (e8, d8)) in
                PointSet::new(self.face_center(face_type).vertex_orbits().into_iter()).iter()
            {
                write!(writer, "{}", vertices.len())?;

                for &vertex in &vertices {
                    write!(writer, " {}", point_sets[0].index(vertex * e8 * d8))?;
                }
                write!(writer, "\n")?;
            }
        }
        write!(writer, "\n")?;

        for i in 3..8 {
            write!(writer, "# {i}-faces\n")?;
            for &face_type in face_types[i].iter() {
                println!("{i}-faces: {face_type:?}");
                let mut subfaces = FxHashSet::from_iter([]);
                for &subface_type in &face_types[i - 1] {
                    if face_type.contains(subface_type) {
                        let mut stack =
                            VecDeque::from_iter([self.face_center(subface_type).vertex()]);
                        while let Some(v) = stack.pop_front() {
                            if subfaces.insert(v) {
                                for mirror in face_type.mirrors() {
                                    stack.push_back(v * mirror.mat());
                                }
                            }
                        }
                    }
                }

                // already did this once, is that ok
                for (_point, (e8, d8)) in
                    PointSet::new(self.face_center(face_type).vertex_orbits().into_iter()).iter()
                {
                    write!(writer, "{}", subfaces.len())?;
                    for &subface in &subfaces {
                        write!(writer, " {}", point_sets[i - 1].index(subface * e8 * d8))?;
                    }
                    write!(writer, "\n")?;
                }
            }
            write!(writer, "\n")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point_set_test(point: [i16; 8]) {
        let point = Point::new(point.into());
        let point_set = PointSet::new([(point, E8::identity())].into_iter());

        assert_eq!(point_set.iter().count() as u64, point_set.len());
        for (i, (point2, (e8, d8))) in point_set.iter().enumerate() {
            assert!(d8.signs_even());
            assert_eq!(point * e8 * d8, point2);
            assert_eq!(point_set.index(point2), i as u64);
        }
    }

    #[test]
    fn point_set_e8_d8_12223344() {
        point_set_test([1, 2, 2, 2, 3, 3, 4, 4]);
    }

    #[test]
    fn point_set_e8_d8_11111111() {
        point_set_test([1, 1, 1, 1, 1, 1, 1, 1]);
    }

    #[test]
    fn point_set_e8_d8_00000022() {
        point_set_test([0, 0, 0, 0, 0, 0, 2, 2]);
    }
}
