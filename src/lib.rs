pub mod backbone;
pub mod external;
pub mod graph;
pub mod io;
pub mod map;
pub mod minimizer;
pub mod molecules;
pub mod overlap;
pub mod protocol;
pub mod repeat;
pub mod report;
pub mod scaffold;

use std::time::Instant;

/// Global timer for consistent elapsed-time logging across the pipeline.
pub struct Timer {
    start: Instant,
}

impl Timer {
    pub fn new() -> Self {
        Self {
            start: Instant::now(),
        }
    }

    pub fn elapsed_secs(&self) -> u64 {
        self.start.elapsed().as_secs()
    }

    pub fn log(&self, msg: &str) {
        log::info!("[{}s] {}", self.elapsed_secs(), msg);
    }
}

impl Default for Timer {
    fn default() -> Self {
        Self::new()
    }
}
