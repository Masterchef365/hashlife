use std::collections::HashMap;
fn main() {

}

type Coord = (i32, i32);
type Unit2x2 = [bool; 4];

/// Macrocell's address in the form [tl, tr, bl, br]
type MacroCellAddress = [usize; 4];

#[derive(Copy, Clone, Debug)]
enum Node {
    Leaf(Unit2x2),
    Branch(MacroCellAddress),
}

struct Engine {
    nodes: Vec<Node>,
    cache: HashMap<[usize; 4], usize>,
}

impl Engine {
    fn query(&mut self, (x, y): Coord, n: usize) -> usize {
        if n == 0 {
            todo!("Query bottom");
        }

        let side_len = 1 << n;
        let tl = self.query((x, y), n - 1);
        let tr = self.query((x, y + side_len), n - 1);
        let bl = self.query((x + side_len, y), n - 1);
        let br = self.query((x + side_len, y + side_len), n - 1);

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
        self.nodes.push(result);
        self.cache.insert(macro_cell, result_idx);
        result_idx
    }
}

fn solve_4x4(tl: Unit2x2, tr: Unit2x2, bl: Unit2x2, br: Unit2x2) -> Unit2x2 {
    todo!()
}
