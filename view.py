import pandas
from chainconsumer import Chain, ChainConsumer, Truth, PlotConfig

filename = "./mock_2/result.pkl"
df = pandas.read_pickle(filename)
columns = ["fs8_eff","Om","MB[1]","alpha", "beta_B",'sigma_int[1]','sigma_int[2]',  "sigma_v"]
df=df[columns]
c = ChainConsumer()
c.add_chain(Chain(samples=df, name="An Example Contour"))  
c.add_truth(Truth(location={"fs8_eff":1,"Om":0.3089,"MB[1]":-19.12,"alpha":0.14, "beta_B":2.9,"sigma_int[1]":0.02,"sigma_int[2]":0.05}))
c.set_plot_config(
    PlotConfig(
        summary_font_size=12, label_font_size=20, legend_kwargs={'prop':{'size':20}}, labels={"fs8_eff":r"$f\sigma_{8,\text{eff}}$","Om": r"$\Omega_M$", "MB[1]":r"$M_B[1]$","alpha":r"$\alpha$","beta_B":r"$\beta$","sigma_int[1]":r"$\sigma_\text{int}[1]$", "sigma_int[2]":r"$\sigma_\text{int}[2]$", "sigma_v":r"$\sigma_v$"}
    )
)
fig = c.plotter.plot(columns=columns)
fig.savefig("temp.png")
