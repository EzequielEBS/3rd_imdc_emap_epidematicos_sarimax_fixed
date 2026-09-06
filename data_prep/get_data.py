import os
import mosqlient

api_key = os.getenv("MOSQLIMATE_API_KEY")

df = mosqlient.get_infodengue(
    api_key=api_key,
    disease="dengue",
    start_date="2022-01-01",
    end_date="2023-01-01",
    uf="AL",
)

print(df.head())