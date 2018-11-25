This repository contains:
- a working draft of a manuscript describing two variations on an IBD detection algorithm commonly used in genetic genealogy applications; the algorithms aim to offer:
  - increased sensitivity and/or reduced false positive rate
  - more accurate segment endpoint detection (with confidence intervals on endpoint location in the case of the second algorithm)
  - comparable algorithmic simplicity and (low) computational expense, compared to the current standard
- a crude python implementation of one of the algorithms, which is simple enough that it could be readily tailored to particular implementations, and could make use of general-purpose GPU hardware
