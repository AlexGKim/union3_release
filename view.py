import pandas
from chainconsumer import Chain, ChainConsumer, Truth

filename = "./output/result.pkl"
df = pandas.read_pickle(filename)
columns = ["fs8_eff","Om","MB[1]","alpha", "beta_B",'sigma_int[1]','sigma_int[2]',  "sigma_v"]
df=df[columns]
c = ChainConsumer()
c.add_chain(Chain(samples=df, name="An Example Contour"))  
c.add_truth(Truth(location={"fs8_eff":1,"Om":0.3089,"MB[1]":-19.053,"alpha":0.14, "beta_B":2.9,"sigma_int[1]":0.02,"sigma_int[2]":0.05}))
fig = c.plotter.plot(columns=columns)
fig.savefig("temp.png")