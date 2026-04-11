import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def plot_enrichment(df:pd.DataFrame,outpng:str, top_n:int=15):
    # 过滤掉 Odds Ratio 为 inf 的项以便绘图，或给它们一个固定最大值
    plot_df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['odds_ratio']).copy()
    plot_df = plot_df.sort_values('odds_ratio', ascending=False).head(top_n)
    
    plt.figure(figsize=(10, 8))
    sns.barplot(
        data=plot_df, 
        x='odds_ratio', 
        y='repeat', 
        palette='viridis'
    )
    
    # 添加一条 OR=1 的参考线
    plt.axvline(x=1, color='red', linestyle='--')
    plt.title(f'Top {top_n} Enriched Repeats in PlaB Group (max_div=3)', fontsize=14)
    plt.xlabel('Odds Ratio (Enrichment)', fontsize=12)
    plt.ylabel('Repeat Subfamily', fontsize=12)
    plt.grid(axis='x', linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.savefig(outpng)

if __name__ == "__main__":
    df = pd.read_csv("/disk5/luosg/Totipotent20251031/a.csv",sep="\t")
    plot_enrichment(df,"a.png")