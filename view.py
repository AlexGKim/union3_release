import pandas
from chainconsumer import Chain, ChainConsumer

filename = "./output/result.pkl"
df = pandas.read_pickle(filename)
columns = ["fs8_eff","Om","MB[1]","alpha", "beta_B",'sigma_int[1]','sigma_int[2]',  "sigma_v"]
df=df[columns]
c = ChainConsumer()
c.add_chain(Chain(samples=df, name="An Example Contour"))
fig = c.plotter.plot(columns=columns)
fig.save("temp.png")