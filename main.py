import pandas as pd
import matplotlib.pyplot as plt

# Read CSV file
# Assume CSV format:
# index,value
# 0,1.2
# 1,3.4
data = pd.read_csv("data/IoT-sensor/PM10Dust.csv")

# Extract columns
x = data.iloc[:5000, 0]
y = data.iloc[:5000, 1]

# Plot
plt.figure(figsize=(8, 4))
plt.plot(x, y)

# Labels
plt.xlabel("Index")
plt.ylabel("Value")
plt.title("CSV Data Plot")

# Save to PNG
plt.savefig("plot.png", dpi=300, bbox_inches="tight")

# Optional: show plot
# plt.show()