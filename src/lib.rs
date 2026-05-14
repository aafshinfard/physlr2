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

/// Global timer for consistent elapsed-time and memory logging across the pipeline.
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

    /// Log a message with elapsed time and current/peak memory usage.
    pub fn log(&self, msg: &str) {
        let (current_mb, peak_mb) = memory_usage_mb();
        log::info!(
            "[{}s, RSS={:.0} MB, peak={:.0} MB] {}",
            self.elapsed_secs(),
            current_mb,
            peak_mb,
            msg
        );
    }
}

impl Default for Timer {
    fn default() -> Self {
        Self::new()
    }
}

/// Read current and peak RSS from /proc/self/status (Linux).
/// Returns (current_rss_mb, peak_rss_mb). Returns (0, 0) on non-Linux.
fn memory_usage_mb() -> (f64, f64) {
    #[cfg(target_os = "linux")]
    {
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            let mut current_kb = 0u64;
            let mut peak_kb = 0u64;
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    current_kb = line
                        .split_whitespace()
                        .nth(1)
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(0);
                } else if line.starts_with("VmHWM:") {
                    peak_kb = line
                        .split_whitespace()
                        .nth(1)
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(0);
                }
            }
            return (current_kb as f64 / 1024.0, peak_kb as f64 / 1024.0);
        }
    }
    (0.0, 0.0)
}
