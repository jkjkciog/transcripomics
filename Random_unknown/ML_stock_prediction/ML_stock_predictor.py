import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

# Load historical stock price data
data = pd.read_csv('historical_stock_prices.csv')
data['Date'] = pd.to_datetime(data['Date'])
data.set_index('Date', inplace=True)

# Feature engineering: create moving averages
data['MA_50'] = data['Close'].rolling(window=50).mean()
data['MA_200'] = data['Close'].rolling(window=200).mean()
data.dropna(inplace=True)

# Normalization
scaler = MinMaxScaler()
scaled_data = scaler.fit_transform(data)

# Splitting data into training and test sets
train_size = int(len(scaled_data) * 0.8)
train, test = scaled_data[:train_size], scaled_data[train_size:]


def stochastic_oscillator(data, lookback):
    data['Lowest_Low'] = data['Low'].rolling(window=lookback).min()
    data['Highest_High'] = data['High'].rolling(window=lookback).max()
    data['%K'] = (data['Close'] - data['Lowest_Low']) / (data['Highest_High'] - data['Lowest_Low']) * 100
    data['%D'] = data['%K'].rolling(window=3).mean()
    return data

# Example usage with a DataFrame `df`
df = stochastic_oscillator(df, lookback=14)