use std::collections::HashMap;
fn main() {
    let data = [
        false, false, true, false,
        true, false, true, false,
        false, true, true, false,
        false, false, false, false,
    ];
    let rect = ((0, 0), (3, 3));
    let mut engine = Engine::default();
    let root = engine.query((0, 0), 2, rect, &data);
    dbg!(engine.nodes[root]);
    //dbg!(engine.cache);
}

type Coord = (i32, i32);
// Inclusive rectangle in the form (min, max). Inclusive means min..=max
type Rect = (Coord, Coord);
type Unit2x2 = [bool; 4];

/// Macrocell's address in the form [tl, tr, bl, br]
type MacroCellAddress = [usize; 4];

#[derive(Copy, Clone, Debug)]
enum Node {
    Leaf(Unit2x2),
    Branch(MacroCellAddress),
}

#[derive(Default)]
struct Engine {
    nodes: Vec<Node>,
    cache: HashMap<[usize; 4], usize>,
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
        let width = x2 - x1 + 1;
        let idx = dx + dy * width;
        input[idx as usize]
    })
}

impl Engine {
    fn query(&mut self, (x, y): Coord, n: usize, input_rect: Rect, input: &[bool]) -> usize {
        //dbg!(n, (x, y));

        if n == 0 {
            panic!("Leaf nodes are 2x2!");
        }

        if n == 1 {
            let bits = [(x, y), (x + 1, y), (x, y + 1), (x + 1, y + 1)]
                .map(|pos| sample_input_rect(pos, input_rect, input).unwrap_or(false));
            //dbg!(bits);
            // TODO: Memoize 2x2 patterns!
            let idx = self.nodes.len();
            self.nodes.push(Node::Leaf(bits));
            return idx;
        }

        let side_len = 1 << n - 1;
        let tl = self.query((x, y), n - 1, input_rect, input);
        let tr = self.query((x + side_len, y), n - 1, input_rect, input);
        let bl = self.query((x, y + side_len), n - 1, input_rect, input);
        let br = self.query((x + side_len, y + side_len), n - 1, input_rect, input);

        self.calc_result([tl, tr, bl, br])
    }

    fn calc_result(&mut self, macro_cell: MacroCellAddress) -> usize {
        if let Some(&result) = self.cache.get(&macro_cell) {
            return result;
        }

        let sub_nodes = macro_cell.map(|idx| self.nodes[idx]);

        let result = match sub_nodes {
            [
                Node::Branch(tl @ [a, b, c, d]), // Top left
                Node::Branch(tr @ [e, f, g, h]), // Top right
                Node::Branch(bl @ [i, j, k, l]), // Bottom left
                Node::Branch(br @ [m, n, o, p]) // Bottom right
            ] =>
            {
                /* 
                    | TL | TR |
                    +----+----|
                    | BL | BR |
                */

                /*
                   | A B | E F |
                   | C D | G G |
                   +-----+-----+
                   | I J | M N |
                   | K L | O P |
                */

                // Top inner row
                let q = self.calc_result(tl);
                let r = self.calc_result([b, e, d, g]);
                let s = self.calc_result(tr);

                // Middle inner row
                let t = self.calc_result([c, d, i, j]);
                let u = self.calc_result([d, g, j, m]);
                let v = self.calc_result([g, h, m, n]);

                // Bottom inner row
                let w = self.calc_result(bl);
                let x = self.calc_result([j, m, l, o]);
                let y = self.calc_result(br);

                /*
                   | Q R S |
                   | T U V |
                   | W X Y |
                */

                // Solve inner 9x9
                let result_tl = self.calc_result([q, r, t, u]);
                let result_tr = self.calc_result([r, s, u, v]);

                let result_bl = self.calc_result([t, u, w, x]);
                let result_br = self.calc_result([u, v, x, y]);

                Node::Branch([result_tl, result_tr, result_bl, result_br])
            }
            [Node::Leaf(tl), Node::Leaf(tr), Node::Leaf(bl), Node::Leaf(br)] => {
                Node::Leaf(solve_4x4(tl, tr, bl, br))
            }
            _ => unreachable!("Mixed cell sizes!"),
        };

        let result_idx = self.nodes.len();
        //dbg!(result);
        self.nodes.push(result);
        self.cache.insert(macro_cell, result_idx);
        result_idx
    }
}

fn solve_4x4(
    [a, b, c, d]: Unit2x2,
    [e, f, g, h]: Unit2x2,
    [i, j, k, l]: Unit2x2,
    [m, n, o, p]: Unit2x2,
) -> Unit2x2 {
    [
        solve_3x3([a, b, e, c, d, g, i, j, m]),
        solve_3x3([b, e, f, d, g, h, j, m, n]),
        solve_3x3([c, d, g, i, j, m, k, l, o]),
        solve_3x3([d, g, h, j, m, n, l, o, p]),
    ]
}

fn solve_3x3([a, b, c, d, e, f, g, h, i]: [bool; 9]) -> bool {
    let count = [a, b, c, d, f, g, h, i].into_iter().filter(|&x| x).count();
    gol_rules(e, count)
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
                    false, false, //.
                    true, false //.
                ],
                [
                    true, false, //.
                    true, false, //.
                ],
                [
                    false, true, //.
                    false, false, //.
                ],
                [
                    true, false, //.
                    false, false, //.
                ]
            ),
            [
                false, true, //.
                true, true, //.
            ]
        );
    }
}
