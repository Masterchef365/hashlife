use std::collections::HashMap;
fn main() {
    let data = [
        false, false, true, false, //.
        true, false, true, false, //.
        false, true, true, false, //.
        false, false, false, false, //.
    ];
    let rect = ((0, 0), (3, 3));
    let mut engine = Engine::new();
    let n = 2;
    let root = engine.query((0, 0), n, rect, &data);

    let cannon = engine.cannon(n, root);
    for row in cannon.chunks_exact(1 << n) {
        for &elem in row {
            print!("{} ", if elem { '#' } else { '_' });
        }
        println!()
    }

    /*
    dbg!(&engine.lookup);
    dbg!(&engine.cache);
    dbg!(root);
    dbg!(engine.lookup[root]);
    */
}

type Coord = (i32, i32);
// Inclusive rectangle in the form (min, max). Inclusive means min..=max
type Rect = (Coord, Coord);

/// Macrocell's address in the form [tl, tr, bl, br]
type MacroCell = [usize; 4];

struct Engine {
    cache: HashMap<[usize; 4], usize>,
    lookup: Vec<[usize; 4]>,
}

fn inside_rect((x, y): Coord, ((x1, y1), (x2, y2)): Rect) -> bool {
    debug_assert!(x1 < x2);
    debug_assert!(y1 < y2);
    x >= x1 && x <= x2 && y >= y1 && y <= y2
}

fn sample_input_rect(
    pos @ (x, y): Coord,
    rect @ ((x1, y1), (x2, _)): Rect,
    input: &[bool],
) -> Option<bool> {
    inside_rect(pos, rect).then(|| {
        let (dx, dy) = (x - x1, y - y1);
        let width = x2 - x1 + 1; // Plus one since the rect is inclusive
        let idx = dx + dy * width;
        input[idx as usize]
    })
}

fn subcoords((x, y): Coord, n: usize) -> [Coord; 4] {
    let side_len = 1 << n;
    [
        (x, y),
        (x + side_len, y),
        (x, y + side_len),
        (x + side_len, y + side_len),
    ]
}

impl Engine {
    pub fn new() -> Self {
        Self {
            lookup: vec![[0; 4], [usize::MAX; 4]],
            cache: HashMap::new(),
        }
    }

    pub fn cannon(&self, n: usize, root: usize) -> Vec<bool> {
        let side_len = 1 << n;
        let mut data = vec![false; side_len * side_len];
        self.cannon_rec((0, 0), n, &mut data, side_len, root);
        data
    }

    fn cannon_rec(&self, corner: Coord, n: usize, buf: &mut [bool], width: usize, node: usize) {
        debug_assert_ne!(n, 0);
        let result = self.lookup[node];

        dbg!(n, corner, result);

        if n == 1 {
            for ((x, y), val) in subcoords(corner, 0).into_iter().zip(result) {
                buf[(x as usize + y as usize * width)] = match val {
                    0 => false,
                    1 => true,
                    other => panic!("N = 1 but {} is not a bit!", other),
                };
            }
        } else {
            for (sub_corner, node) in subcoords(corner, n - 1).into_iter().zip(result) {
                self.cannon_rec(sub_corner, n - 1, buf, width, node);
            }
        }
    }

    fn query(&mut self, corner: Coord, n: usize, input_rect: Rect, input: &[bool]) -> usize {
        if n == 0 {
            panic!("Leaf nodes are 2x2!");
        }

        if n == 1 {
            let bits = subcoords(corner, 0)
                .map(|pos| sample_input_rect(pos, input_rect, input).unwrap_or(false) as usize);

            dbg!(n, bits);
            return match self.cache.get(&bits) {
                Some(&idx) => idx,
                None => dbg!(self.push_cell(bits)),
            };
        }

        let sub_squares = subcoords(corner, n - 1)
            .map(|sub_corner| self.query(sub_corner, n - 1, input_rect, input));

        self.calc_result(sub_squares, n)
    }

    fn calc_result(&mut self, macro_cell: MacroCell, level: usize) -> usize {
        eprintln!("CALC {}", level);
        if let Some(&result) = self.cache.get(&macro_cell) {
            eprintln!("CACHED {:?} -> {}", macro_cell, result);
            return result;
        }

        let [
            tl @ [_, b, c, d], // Top left
            tr @ [e, _, g, h], // Top right
            bl @ [i, j, _, l], // Bottom left
            br @ [m, n, o, _] // Bottom right
        ] = macro_cell.map(|idx| self.lookup[idx]);

        if level == 1 {
            return self.push_cell(macro_cell);
        }

        if level == 2 {
            // Check if we've encountered this node before
            let soln = solve_4x4(tl, tr, bl, br);
            eprintln!("PRIMITIVE {:?} -> {:?}", macro_cell, soln);
            match self.cache.get(&soln) {
                Some(&idx) => idx,
                None => dbg!(self.push_cell(soln)),
            }
        } else {
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
            let q = self.calc_result(tl, level - 2);
            let r = self.calc_result([b, e, d, g], level - 2);
            let s = self.calc_result(tr, level - 2);

            // Middle inner row
            let t = self.calc_result([c, d, i, j], level - 2);
            let u = self.calc_result([d, g, j, m], level - 2);
            let v = self.calc_result([g, h, m, n], level - 2);

            // Bottom inner row
            let w = self.calc_result(bl, level - 2);
            let x = self.calc_result([j, m, l, o], level - 2);
            let y = self.calc_result(br, level - 2);

            /*
            | Q R S |
            | T U V |
            | W X Y |
            */

            // Solve inner 9x9
            let result_tl = self.calc_result([q, r, t, u], level - 2);
            let result_tr = self.calc_result([r, s, u, v], level - 2);

            let result_bl = self.calc_result([t, u, w, x], level - 2);
            let result_br = self.calc_result([u, v, x, y], level - 2);

            let result = [result_tl, result_tr, result_bl, result_br];

            eprintln!("CONSTRUCT {:?} -> {:?}", macro_cell, result);
            dbg!(self.push_cell(result))
        }
    }

    fn push_cell(&mut self, cell: MacroCell) -> usize {
        let result_idx = self.lookup.len();
        self.lookup.push(cell);
        self.cache.insert(cell, result_idx);
        result_idx
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
}
