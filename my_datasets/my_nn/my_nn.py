import torch
import torch.nn as nn
import torch.nn.functional as F
import json
from torch.utils.data import Dataset, DataLoader

class my_dataset(Dataset):

    def __init__(self, root, train=True):
        super().__init__()
        self.train = train
        data = []
        
        if self.train:
            data_path = root + 'train.txt'
        else:
            data_path = root + 'test.txt'
            
        with open(data_path, 'r') as f:
            for line in f:
                item = json.loads(line)
                data.append(item)

        self.sequences = [x['seq'] for x in data]
        self.labels = [x['tag'] for x in data]

    def __getitem__(self, index):
        seq = self.sequences[index]
        label = self.labels[index]
        return seq, label

    def __len__(self):
        return len(self.sequences)


class CNN_net(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Conv1d(1, 1, 6)
        self.pool = nn.AvgPool1d(2)
        self.fc = nn.Linear(5 * 5 * 2, 2)

    def forward(self, x):
        out = self.conv(x)
        out = F.relu(out)
        out = self.pool(out)
        out = out.view(-1, self.num_flat_features(out))
        out = self.fc(out)
        return out

    def num_flat_features(self, x):
        size = x.size()[1:]
        num_features = 1
        for s in size:
            num_features *= s
        return num_features


train_dataset = my_dataset('Bioseq/data/Prom', train=True)
test_dataset = my_dataset('Bioseq/data/Prom', train=False)

train_loader = DataLoader(dataset=train_dataset, batch_size=64, shuffle=True)
test_loader = DataLoader(dataset=test_dataset, batch_size=64, shuffle=False)

model = CNN_net()
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=0.01)

# Training loop
for epoch in range(2):  # loop over the dataset multiple times
    for i, data in enumerate(train_loader, 0):
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

print('Finished Training')
