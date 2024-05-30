SELECT DISTINCT 
	[NHSNumber],
	[FinancialYear],
	[AdmissionDate],
	[DischargeDate],
	[AgeOnAdmission],
	[DiagnosisCode],
	[DiagnosisDescription],
	--[OSLAUA] AS LocalAuthCode,
	[Cost]
  FROM 
	[EAT_Reporting_BSOL].[SUS].[VwInpatientEpisodesPatientGeography] AS A
  LEFT JOIN 
	[EAT_Reporting_BSOL].[SUS].[VwInpatientEpisodesDiagnosisRelational] AS B
  ON 
	A.[EpisodeId] = B.[EpisodeId]
  WHERE 
    ([DiagnosisCode] LIKE 'A80%' OR
	 [DiagnosisCode] LIKE 'G14%' OR
	 [DiagnosisCode] LIKE 'B91%'
	) AND
	([AdmissionDate] >= '2013-04-01' AND [AdmissionDate] < '2024-04-01') AND
	Cost IS NOT NULL 


SELECT DISTINCT 
	[NHSNumber],
	[AdmissionDate],
	[OSLAUA]
  FROM 
	[EAT_Reporting_BSOL].[SUS].[VwInpatientEpisodesPatientGeography] AS A
  LEFT JOIN 
	[EAT_Reporting_BSOL].[SUS].[VwInpatientEpisodesDiagnosisRelational] AS B
  ON 
	A.[EpisodeId] = B.[EpisodeId]
  WHERE 
    ([DiagnosisCode] LIKE 'A80%') AND
	([AdmissionDate] >= '2020-04-01' AND [AdmissionDate] < '2024-04-01')