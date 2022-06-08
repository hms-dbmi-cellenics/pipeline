{{/* Generate a template that can be used for both assigned and unassigned xperiments */}}
# go templating language
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
            value: "{{ .Values.myAccount.accountId }}"
          - name: AWS_DEFAULT_REGION
            value: "{{ .Values.myAccount.region }}"
        volumeMounts:
        - name: podinfo
          mountPath: /etc/podinfo
        resources:
          requests:
            memory: "{{ .Values.memoryRequest }}"
      volumes:
      - name: podinfo
        downwardAPI:
          items:
            - path: "labels"
              fieldRef:
                fieldPath: metadata.labels
{{- end -}}
