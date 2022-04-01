//! Naive wasm-friendly fastq parser
//!
//! Doesn't allocate new buffers for each record

struct FastqIter<R> {
    src: BufReader<R>,
    name: Rc<String>,
    seq: Rc<String>,
    sep: Rc<String>,
    qual: Rc<String>,
}

impl<R: Read> FastqIter<R> {
    fn new(src: BufReader<R>) -> FastqIter<R> {
        FastqIter {
            src: src,
            name: Rc::new(String::new()),
            seq: Rc::new(String::new()),
            sep: Rc::new(String::new()),
            qual: Rc::new(String::new()),
        }
    }
}

struct Record {
    seq: Rc<String>,
}

impl<R: Read> Iterator for FastqIter<R> {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        let name = match Rc::get_mut(&mut self.name) {
            Some(name) => {
                name.clear();
                name
            }
            None => {
                self.name = Rc::new(String::new());
                Rc::make_mut(&mut self.name)
            }
        };
        let seq = match Rc::get_mut(&mut self.seq) {
            Some(seq) => {
                seq.clear();
                seq
            }
            None => {
                self.seq = Rc::new(String::new());
                Rc::make_mut(&mut self.seq)
            }
        };
        let sep = match Rc::get_mut(&mut self.sep) {
            Some(sep) => {
                sep.clear();
                sep
            }
            None => {
                self.sep = Rc::new(String::new());
                Rc::make_mut(&mut self.sep)
            }
        };
        let qual = match Rc::get_mut(&mut self.qual) {
            Some(qual) => {
                qual.clear();
                qual
            }
            None => {
                self.qual = Rc::new(String::new());
                Rc::make_mut(&mut self.qual)
            }
        };

        self.src.read_line(name);
        self.src.read_line(seq);
        self.src.read_line(sep);
        match self.src.read_line(qual) {
            Ok(0) => None,
            Ok(_) => Some(Record {
                seq: Rc::clone(&self.seq),
            }),
            Err(_) => None,
        }
    }
}

