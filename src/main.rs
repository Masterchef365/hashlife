use ahash::{AHasher, RandomState};
use std::collections::HashMap;
mod io;
use io::load_rle;

fn main() {
    let mut args = std::env::args().skip(1);
    let n = args.next().unwrap_or("2".into()).parse().unwrap();
    let rle_path = args.next();

    let (data, width) = if let Some(rle_path) = rle_path {
        load_rle(rle_path).expect("IO Error")
    } else {
        (
            vec![
                false, false, true, false, //.
                true, false, true, false, //.
                false, true, true, false, //.
                false, false, false, false, //.
            ],
            4,
        )
    };

    let height = data.len() / width;

    let half_n_sq = 1i64 << n - 1;
    let pos @ (x, y) = (half_n_sq, half_n_sq);
    let rect = (pos, (x + width as i64, y + height as i64));

    let mut engine = Engine::new();

    let qt = std::time::Instant::now();
    let root = engine.query((0, 0), n, rect, &data);
    let qt = qt.elapsed();

    let root = engine.lookup[root].1.unwrap();

    /*for (idx, elem) in engine.lookup.iter().enumerate() {
        println!("{}, {:?}", idx, elem);
    }
    dbg!(root);*/

    println!("Time: {}ms", qt.as_secs_f32() * 1000.);

    let out_width = 300;
    let half = 1i64 << n - 2;
    let rect = extend_rect(((half, half), (half, half)), out_width);

    let cannon = engine.cannon(n - 1, root, rect);
    /*
    for row in cannon.chunks_exact(1 << n - 1) {
        for &elem in row {
            print!("{} ", if elem { '#' } else { '_' });
        }
        println!()
    }
    */

    let (width, _) = rect_dimensions(rect);
    let image_data: Vec<u8> = cannon.into_iter().map(|b| [b as u8 * 255; 3]).flatten().collect();
    io::write_ppm("out.ppm", &image_data, width as _).expect("Image IO error");
}

type Coord = (i64, i64);
// Inclusive rectangle in the form (min, max). Inclusive means min..=max
type Rect = (Coord, Coord);

/// Macrocell's address in the form [tl, tr, bl, br]
type MacroCell = [usize; 4];

struct Engine {
    /// Maps a set of quadrants to the parent macrocell
    cache: HashMap<[usize; 4], usize, RandomState>,
    /// Macrocells, each lists it's quadrants and it's result cell (if present)
    lookup: Vec<([usize; 4], Option<usize>)>,
}

/// Check if the given point is inside the given rectangle
fn inside_rect((x, y): Coord, ((x1, y1), (x2, y2)): Rect) -> bool {
    debug_assert!(x1 < x2);
    debug_assert!(y1 < y2);
    x >= x1 && x < x2 && y >= y1 && y < y2
}

/// Sample at `pos` from the given `input` buffer positioned at `rect`
fn sample_rect(
    pos @ (x, y): Coord,
    rect @ ((x1, y1), (x2, _)): Rect,
) -> Option<usize> {
    // debug_assert_eq!(input.len(), (x2 - x1) * (y2 - y1)) // TODO:
    inside_rect(pos, rect).then(|| {
        let (dx, dy) = (x - x1, y - y1);
        let width = x2 - x1; // Plus one since the rect is inclusive
        (dx + dy * width) as usize
    })
}

/// Given a corner position at `(x, y)`, and a size `n` return the corner positions of the four
/// quadrants that make up the macrocell.
fn subcoords((x, y): Coord, n: usize) -> [Coord; 4] {
    let side_len = 1 << n;
    [
        (x, y),
        (x + side_len, y),
        (x, y + side_len),
        (x + side_len, y + side_len),
    ]
}

/// Returns the (width, height) of the given rect
fn rect_dimensions(((x1, y1), (x2, y2)): Rect) -> (i64, i64) {
    (x2 - x1, y2 - y1)
}

impl Engine {
    pub fn new() -> Self {
        Self {
            lookup: vec![
                ([0; 4], Some(0)),       // The zero cell always results in a zero cell
                ([usize::MAX; 4], None), // The one cell does not map to anything!
            ],
            cache: [([0; 4], 0)].into_iter().collect(),
        }
    }

    pub fn cannon(&self, n: usize, root: usize, rect: Rect) -> Vec<bool> {
        let (width, height) = rect_dimensions(rect);
        let mut data = vec![false; (width * height) as usize];
        self.cannon_rec((0, 0), n, &mut data, rect, root);
        data
    }

    fn cannon_rec(&self, corner: Coord, n: usize, buf: &mut [bool], rect: Rect, cell: usize) {
        debug_assert_ne!(n, 0);
        let (quadrants, _) = self.lookup[cell];

        if zero_input(corner, n, rect) {
            return;
        }

        if n == 1 {
            for (pos, val) in subcoords(corner, 0).into_iter().zip(quadrants) {
                if let Some(idx) = sample_rect(pos, rect) {
                    buf[idx] = match val {
                        0 => false,
                        1 => true,
                        other => panic!("N = 1 but {} is not a bit!", other),
                    };
                }
            }
        } else {
            for (sub_corner, node) in subcoords(corner, n - 1).into_iter().zip(quadrants) {
                self.cannon_rec(sub_corner, n - 1, buf, rect, node);
            }
        }
    }

    /// Calculate the result square of the given area of input with time stamp
    fn query(&mut self, corner: Coord, n: usize, input_rect: Rect, input: &[bool]) -> usize {
        // Short circuit for zeroes
        if zero_input(corner, n, input_rect) {
            return 0;
        }

        // Return the input pixel at the given coordinates
        if n == 0 {
            return sample_rect(corner, input_rect)
                .map(|idx| input[idx])
                .unwrap_or(false) as usize;
        }

        // Calculate which macrocell we are in
        let macro_cell = subcoords(corner, n - 1)
            .map(|sub_corner| self.query(sub_corner, n - 1, input_rect, input));

        // Calculate the result for this macro cell
        let result_idx = (n >= 2).then(|| self.calc_result(macro_cell, n));

        // Push the new macrocell
        self.push_cell(macro_cell, result_idx)
    }

    /// Calculate the result square of the given macro cell
    fn calc_result(&mut self, macro_cell: MacroCell, level: usize) -> usize {
        if level <= 1 {
            panic!("We can't compute results for 1x1s or 2x2s!");
        }

        // Check if we already know the result
        if let Some(result) = self
            .cache
            .get(&macro_cell)
            .and_then(|&idx| self.lookup[idx].1)
        {
            return result;
        }

        // Deconstruct the quadrants of the macrocell
        let [
            (tl @ [_, b, c, d], _), // Top left
            (tr @ [e, _, g, h], _), // Top right
            (bl @ [i, j, _, l], _), // Bottom left
            (br @ [m, n, o, _], _) // Bottom right
        ] = macro_cell.map(|idx| self.lookup[idx]);

        let result = if level == 2 {
            // Compute the 4x4 from scratch. If the resulting 2x2 is already in cache, return that
            solve_4x4(tl, tr, bl, br)
        } else {
            // We need to compute the result from composite 3x3 and 2x2
            /*
            | TL | TR |
            +----+----+
            | BL | BR |
            */

            /*
            | _ B | E _ |
            | C D | G H |
            +-----+-----+
            | I J | M N |
            | _ L | O _ |
            */

            // Top inner row
            let q = self.calc_result(tl, level - 1);
            let r = self.calc_result([b, e, d, g], level - 1);
            let s = self.calc_result(tr, level - 1);

            // Middle inner row
            let t = self.calc_result([c, d, i, j], level - 1);
            let u = self.calc_result([d, g, j, m], level - 1);
            let v = self.calc_result([g, h, m, n], level - 1);

            // Bottom inner row
            let w = self.calc_result(bl, level - 1);
            let x = self.calc_result([j, m, l, o], level - 1);
            let y = self.calc_result(br, level - 1);

            /*
            | Q R S |
            | T U V |
            | W X Y |
            */

            // Solve inner 3x3
            let result_tl = self.calc_result([q, r, t, u], level - 1);
            let result_tr = self.calc_result([r, s, u, v], level - 1);

            let result_bl = self.calc_result([t, u, w, x], level - 1);
            let result_br = self.calc_result([u, v, x, y], level - 1);

            [result_tl, result_tr, result_bl, result_br]
        };

        let result_idx = self.push_cell(result, None);

        self.push_cell(macro_cell, Some(result_idx));

        result_idx
    }

    /// Push a cell and it's result, adding to both lookup and cache
    fn push_cell(&mut self, cell: MacroCell, result: Option<usize>) -> usize {
        match self.cache.get(&cell) {
            Some(&idx) => {
                if let Some(result) = result {
                    let _ = self.lookup[idx].1.insert(result);
                }
                idx
            }
            None => {
                let idx = self.lookup.len();
                self.lookup.push((cell, result));
                self.cache.insert(cell, idx);
                idx
            }
        }
    }
}

fn solve_4x4(
    [a, b, c, d]: MacroCell,
    [e, f, g, h]: MacroCell,
    [i, j, k, l]: MacroCell,
    [m, n, o, p]: MacroCell,
) -> MacroCell {
    [
        solve_3x3([a, b, e, c, d, g, i, j, m]),
        solve_3x3([b, e, f, d, g, h, j, m, n]),
        solve_3x3([c, d, g, i, j, m, k, l, o]),
        solve_3x3([d, g, h, j, m, n, l, o, p]),
    ]
}

fn solve_3x3([a, b, c, d, e, f, g, h, i]: [usize; 9]) -> usize {
    let count = [a, b, c, d, f, g, h, i].into_iter().sum();
    gol_rules(e != 0, count) as usize
}

fn gol_rules(center: bool, neighbors: usize) -> bool {
    match (center, neighbors) {
        (true, n) if (n == 2 || n == 3) => true,
        (false, n) if (n == 3) => true,
        _ => false,
    }
}

fn rect_intersect(a: Rect, b: Rect) -> bool {
    let ((x1a, y1a), (x2a, y2a)) = a;
    let ((x1b, y1b), (x2b, y2b)) = b;
    x1a < x2b && x1b < x2a && y1a < y2b && y1b < y2a
}

fn extend_rect(((x1, y1), (x2, y2)): Rect, w: i64) -> Rect {
    ((x1 - w, y1 - w), (x2 + w, y2 + w))
}

/// Calculates whether or not this square can be anything other than zero, given the input rect
fn zero_input(coord: Coord, n: usize, input_rect: Rect) -> bool {
    if n == 0 {
        return false;
    }

    let time = 1 << n - 1;
    let input_rect = extend_rect(input_rect, time);

    let width = 1 << n;
    let (x, y) = coord;
    !rect_intersect((coord, (x + width, y + width)), input_rect)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_4x4() {
        assert_eq!(
            solve_4x4(
                [
                    0, 0, //.
                    1, 0 //.
                ],
                [
                    1, 0, //.
                    1, 0, //.
                ],
                [
                    0, 1, //.
                    0, 0, //.
                ],
                [
                    1, 0, //.
                    0, 0, //.
                ]
            ),
            [
                0, 1, //.
                1, 1, //.
            ]
        );
    }

    #[test]
    fn test_rect_intersect() {
        #[track_caller]
        fn reflexive(a: Rect, b: Rect, expect: bool) {
            assert_eq!(rect_intersect(a, b), expect);
            assert_eq!(rect_intersect(b, a), expect);
        }
        reflexive(((-10, -10), (10, 10)), ((-5, -5), (5, 5)), true);
        reflexive(((-10, -10), (10, 10)), ((0, 0), (5, 5)), true);
        reflexive(((3, 3), (10, 10)), ((0, 0), (5, 5)), true);
        reflexive(((7, 7), (10, 10)), ((0, 0), (5, 5)), false);
        reflexive(((7, 7), (10, 10)), ((0, 0), (5, 5)), false);
        reflexive(((7, 7), (10, 10)), ((-5, -5), (0, 0)), false);
        reflexive(((7, 7), (10, 10)), ((-5, -5), (5, 50)), false);
        reflexive(((0, 0), (10, 10)), ((-5, -5), (5, 50)), true);
        reflexive(((0, 0), (10, 10)), ((-5, -5), (50, 5)), true);

        reflexive(((0, 0), (10, 10)), ((-5, -5), (-2, 50)), false);
        reflexive(((0, 0), (10, 10)), ((-5, -5), (50, -2)), false);
    }
}
