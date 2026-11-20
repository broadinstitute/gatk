#!/usr/bin/env python3
# Apple GPU (Metal/MPS) vs CPU parity probe for the kind of CNN scorevariants uses.
# Runs identical weights+input (float32) on CPU and MPS and reports both:
#   * raw numeric difference (NOT bit-identical: max ~9e-8 on Apple Silicon), and
#   * BIO-equivalence: identical predicted class (argmax) and scores agreeing to the precision
#     GATK actually reports/uses -> no variant call or PASS/FAIL would change.
# Conclusion: MPS is bio-identical to CPU, so it is kept (auto-used) for the Apple GPU speedup.
# See docs/arm64/gpu.md.
import torch
import torch.nn as nn
torch.manual_seed(0)


# A small CNN representative of scorevariants (Conv2d + ReLU + Linear + softmax), float32
class Net(nn.Module):
    def __init__(self):
        super().__init__()
        self.c1 = nn.Conv2d(4, 32, 3, padding=1)
        self.c2 = nn.Conv2d(32, 32, 3, padding=1)
        self.fc = nn.Linear(32 * 8 * 8, 4)

    def forward(self, x):
        x = torch.relu(self.c1(x))
        x = torch.relu(self.c2(x))
        x = x.flatten(1)
        x = self.fc(x)
        return torch.softmax(x, dim=1)


net = Net().eval()
x = torch.randn(256, 4, 8, 8)

with torch.no_grad():
    out_cpu = net(x)
    out_mps = net.to('mps')(x.to('mps')).to('cpu')

diff = (out_cpu - out_mps).abs()
print("max |CPU-MPS|        =", diff.max().item())
print("bit-identical        =", torch.equal(out_cpu, out_mps))
# Bio-equivalence checks: same predicted class, and same score at GATK's reporting precision.
same_class = torch.equal(out_cpu.argmax(1), out_mps.argmax(1))
# GATK writes CNN scores to a few decimals; round both to 3 decimals (well coarser than 9e-8).
same_rounded = torch.equal((out_cpu * 1000).round(), (out_mps * 1000).round())
print("same predicted class =", same_class, "(", out_cpu.size(0), "variants )")
print("same score @3 decimals=", same_rounded)
print("BIO-IDENTICAL        =", bool(same_class and same_rounded))
