import numpy as np
import colour

# --- Step 1: Load cone fundamentals and extract wavelength / LMS arrays ---

cmfs = colour.colorimetry.MSDS_CMFS_LMS["Stockman & Sharpe 2 Degree Cone Fundamentals"]

wavelengths = cmfs.wavelengths          # shape (441,), 390..830 nm @ 1nm
LMS = cmfs.values                       # shape (441, 3), columns = L, M, S

print(f"Wavelengths: {wavelengths[0]:.0f} – {wavelengths[-1]:.0f} nm, "
      f"{len(wavelengths)} samples")
print(f"LMS array shape: {LMS.shape}")

# Peak wavelengths (sanity check: L~560, M~530, S~420-ish)
for i, name in enumerate(["L", "M", "S"]):
    peak = wavelengths[np.argmax(LMS[:, i])]
    print(f"  {name} cone peak: {peak:.0f} nm")

# Save for next steps
np.save("wavelengths.npy", wavelengths)
np.save("LMS.npy", LMS)
print("\nSaved wavelengths.npy and LMS.npy")
