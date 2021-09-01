{{/* Generate a template that can be used for both assigned and unassigned xperiments */}}
{{- define "pipeline.pod-template" -}}
    metadata:
      name: "{{ .Release.Name }}-server"
      labels:
        type: 'pipeline'
        sandboxId: "{{ .Values.sandboxId }}"
    spec:
      restartPolicy: Always
      serviceAccountName: 'deployment-runner'
      containers:
      - name: "{{ .Release.Name }}"
        image: "{{ .Values.pipelineRunner.image }}"
        env:
          - name: CLUSTER_ENV
            value: "{{ .Values.clusterEnv }}"
          - name: SANDBOX_ID
            value: "{{ .Values.sandboxId }}"
          - name: AWS_ACCOUNT_ID
            value: "{{ .Values.awsAccountId }}"
          - name: AWS_DEFAULT_REGION
            value: "{{ .Values.awsRegion }}"
        volumeMounts:
        - name: podinfo
          mountPath: /etc/podinfo
        resources:
          requests:
            memory: "29Gi"
      volumes:
      - name: podinfo
        downwardAPI:
          items:
            - path: "labels"
              fieldRef:
                fieldPath: metadata.labels
{{- end -}}