# Setup

  * Install .NET 6.0 from [here](https://dotnet.microsoft.com/en-us/download/dotnet). Note `CromwellOnAzure` really wants
  .NET version 6.0 and will not run with .NET 7.0.
  * Clone the `CromwellOnAzure` git repo from [here](https://github.com/microsoft/CromwellOnAzure).
  * Build `CromwellOnAzure` with 
```
cd CromwellOnAzure
dotnet build
```

  * Create your deployment of `CromwellOnAzure` in the Variants Azure subscription with:
```
dotnet src/deploy-cromwell-on-azure/bin/Debug/net6.0/deploy-cromwell-on-azure.dll --SubscriptionId 03a6af79-486d-4ada-b0ed-40083610727d --RegionName eastus --MainIdentifierPrefix coa
```
